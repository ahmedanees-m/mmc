"""Step 1: run every structure source through one identical evaluation (PREREG_v4 section 2).

Each source supplies an edge set and nothing else. The edge set is compiled into the
same structural backend, fit with the same optimizer budget, and scored on the same
leave-one-perturbation-out folds, so the only thing that differs between rows of the
output table is where the edges came from.

The oracle arm (S5) and the random null (S3) dominate the cost, because both score
many structures rather than one. Both are dispatched across a process pool; the pool
size is a command-line argument and defaults well below the core count, since the
host runs other work.

    python scripts/step1_comparator.py --module Cytokine_production \
        --condition Stim8hr --sources textbook,random,oracle,linear,mean,zero \
        --workers 12 --seeds 10 --out results/step1_cytokine.json

Requires MMC_ZHU_STORE. The claude source additionally requires a frozen structure
file, and never calls the API from here; S1 structures are frozen upstream so this
script stays deterministic and free.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np

# One thread per worker; the pool owns the parallelism.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false")

from mmc.data import module_data, module_extract, textbook  # noqa: E402
from mmc.eval import compare, identifiability, oracle_search, random_null  # noqa: E402
from mmc.grammar.model_spec import ModelSpec  # noqa: E402

SEARCH_FOLDS = 5          # PREREG_v4 amendment A1
SEARCH_STARTS = 1
SEARCH_MAX_ITER = 80
SCREEN_FOLDS = 2
SCREEN_KEEP = 20
GREEDY_POOL = 120
MAX_EDGES = 30
ANNEAL_STEPS = 800
ANNEAL_BATCH = 24
NESTED_OUTER_FRACTION = 0.30

_W: dict = {}


# ------------------------------ worker side ------------------------------
def _init_worker(genes, perts, observed, de_mask):
    _W["genes"] = list(genes)
    _W["perts"] = list(perts)
    _W["observed"] = np.asarray(observed, float)
    _W["de_mask"] = np.asarray(de_mask, bool)


def _fold_assignment(n: int, n_folds: int, seed: int) -> list[list[int]]:
    if n_folds >= n:
        return [[i] for i in range(n)]
    order = np.random.default_rng(seed).permutation(n)
    return [list(map(int, part)) for part in np.array_split(order, n_folds)]


def _score_edge_set(task):
    """Score one edge set. Returns (mean held-out DE-overlap, mean training loss)."""
    edges, n_folds, n_starts, max_iter, seed, subset = task
    from mmc.compile import structural
    from mmc.eval.holdout import de_overlap
    from mmc.fit import fit_structural

    genes, perts = _W["genes"], _W["perts"]
    obs, de = _W["observed"], _W["de_mask"]
    rows = list(subset) if subset is not None else list(range(len(perts)))
    if not edges or len(rows) < 2:
        return (float("nan"), float("inf"))

    spec = oracle_search.spec_from_edges(genes, edges)
    folds = _fold_assignment(len(rows), n_folds, seed)
    scores, losses = [], []
    for fold in folds:
        held = [rows[i] for i in fold]
        train = [r for r in rows if r not in set(held)]
        if not train:
            continue
        observed = {perts[i]: {genes[j]: float(obs[i, j]) for j in range(len(genes))
                               if genes[j] != perts[i]} for i in train}
        try:
            fits = fit_structural.multi_fit(spec, observed, n_starts=n_starts,
                                            max_iter=max_iter)
        except Exception:
            return (float("nan"), float("inf"))
        losses.append(float(fits[0]["loss"]))
        for i in held:
            pred = np.asarray(structural.knockdown(spec, fits[0]["params"], perts[i]))
            scores.append(de_overlap(pred, obs[i], de[i]))
    if not scores:
        return (float("nan"), float("inf"))
    with np.errstate(invalid="ignore"):
        mean_score = float(np.nanmean(scores)) if not np.all(np.isnan(scores)) else float("nan")
    return (mean_score, float(np.mean(losses)) if losses else float("inf"))


# ------------------------------ driver side ------------------------------
def _make_scorer(pool, n_folds, n_starts, max_iter, seed, subset=None, memo=None):
    """Dispatch a batch of edge sets to the pool, reusing anything already scored.

    The backward pass re-scores structures the forward pass has already seen, and
    annealing revisits states it has left, so a memo on the edge set removes a
    substantial fraction of the work at no cost to the result.
    """
    memo = {} if memo is None else memo

    def score_many(edge_sets):
        keys = [tuple(sorted(es)) for es in edge_sets]
        todo = [k for k in dict.fromkeys(keys) if k not in memo]
        if todo:
            tasks = [(k, n_folds, n_starts, max_iter, seed, subset) for k in todo]
            for k, val in zip(todo, pool.map(_score_edge_set, tasks, chunksize=1)):
                memo[k] = val
        return [memo[k] for k in keys]

    return score_many


def run_oracle(pool, mod, seed, *, subset=None, log=print):
    """One seed of the S5 search. Returns the selected edge set and its trace."""
    trace = oracle_search.SearchTrace()
    score = _make_scorer(pool, SEARCH_FOLDS, SEARCH_STARTS, SEARCH_MAX_ITER, seed, subset,
                         memo={})
    screen = _make_scorer(pool, SCREEN_FOLDS, SEARCH_STARTS, SEARCH_MAX_ITER, seed, subset,
                          memo={})
    pool_pairs = oracle_search.candidate_pool(mod, size=GREEDY_POOL)

    t0 = time.time()
    greedy = oracle_search.greedy_forward_backward(
        score, pool_pairs, max_edges=MAX_EDGES, trace=trace,
        screen_many=screen, screen_keep=SCREEN_KEEP)
    greedy_score = score([frozenset(greedy)])[0][0]
    log(f"    greedy: {len(greedy)} edges, search score {greedy_score:.4f}, "
        f"{time.time() - t0:.0f}s, {len(trace.holdout)} structures scored")

    t1 = time.time()
    annealed = oracle_search.simulated_annealing(
        score, greedy, oracle_search.all_pairs(mod), n_steps=ANNEAL_STEPS,
        batch=ANNEAL_BATCH, max_edges=MAX_EDGES, seed=seed, trace=trace)
    annealed_score = score([frozenset(annealed)])[0][0]
    log(f"    anneal: {len(annealed)} edges, search score {annealed_score:.4f}, "
        f"{time.time() - t1:.0f}s")

    best = annealed if annealed_score >= greedy_score else greedy
    return {
        "edges": sorted(best),
        "greedy_edges": sorted(greedy),
        "search_score_greedy": greedy_score,
        "search_score_annealed": annealed_score,
        "anneal_gain": annealed_score - greedy_score,
        "n_structures_scored": len(trace.holdout),
        "seconds": round(time.time() - t0, 1),
    }, trace


def evaluate_spec(mod, spec, name, seed=0):
    fn = compare.bind_structural(spec, mod)
    return compare.evaluate_source(mod, name, fn, n_edges=len(spec.edges))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--sources", default="textbook,linear,mean,zero,oracle,random")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--n-random", type=int, default=1000)
    ap.add_argument("--random-edges", type=int, default=0,
                    help="edge count for S3; 0 matches the oracle's selected count")
    ap.add_argument("--frozen", default="", help="JSON holding the frozen S1 structure")
    ap.add_argument("--module-def", default="", help="JSON registering a dynamic module")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(args.module, d["regulators"], d["targets"])

    sources = [s.strip() for s in args.sources.split(",") if s.strip()]
    mod = module_data.build_module_data(args.module, args.condition)
    summary = module_data.module_summary(mod)
    print(f"=== Step 1 comparator: {args.module} / {args.condition} ===")
    print(f"{summary['n_genes']} genes, {summary['n_perts']} perturbations, "
          f"{summary['n_de_entries']} DE entries, "
          f"{summary['n_folds_with_de']} folds with a DE gene")

    out: dict = {
        "module": args.module,
        "condition": args.condition,
        "module_summary": summary,
        "protocol": {
            "fit_starts": compare.FIT_STARTS,
            "fit_max_iter": compare.FIT_MAX_ITER,
            "n_boot": compare.N_BOOT,
            "search_folds": SEARCH_FOLDS,
            "search_starts": SEARCH_STARTS,
            "search_max_iter": SEARCH_MAX_ITER,
            "screen_folds": SCREEN_FOLDS,
            "screen_keep": SCREEN_KEEP,
            "greedy_pool": GREEDY_POOL,
            "max_edges": MAX_EDGES,
            "anneal_steps": ANNEAL_STEPS,
            "seeds": args.seeds,
            "n_random": args.n_random,
        },
        "diagnostics": identifiability.diagnostics(mod),
        "sources": {},
        "specs": {},
    }
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    def flush():
        with open(outpath, "w") as f:
            json.dump(out, f, indent=2, default=float)

    results: dict[str, compare.SourceResult] = {}

    # cheap, deterministic sources first so a long run has something on disk early
    if "linear" in sources:
        results["linear"] = compare.evaluate_source(mod, "linear", compare.bind_linear(mod))
        print("  linear done")
    if "mean" in sources:
        results["mean"] = compare.evaluate_source(mod, "mean", compare.bind_mean(mod))
    if "zero" in sources:
        results["zero"] = compare.evaluate_source(mod, "zero", compare.bind_zero(mod))
    flush()

    if "textbook" in sources:
        try:
            spec = textbook.textbook_spec(args.module, mod.genes)
            results["textbook"] = evaluate_spec(mod, spec, "textbook")
            out["specs"]["textbook"] = json.loads(spec.to_json())
            out["textbook_coverage"] = textbook.coverage(args.module, mod.genes)
            print(f"  textbook done ({len(spec.edges)} edges)")
        except KeyError as e:
            print(f"  textbook skipped: {e}")
        flush()

    if "claude" in sources and args.frozen:
        with open(args.frozen) as f:
            frozen = json.load(f)
        spec = ModelSpec.model_validate(frozen["best"]["spec"])
        spec = ModelSpec(genes=list(mod.genes),
                         edges=[e for e in spec.edges
                                if e.regulator in set(mod.genes) and e.target in set(mod.genes)],
                         rules={t: r for t, r in spec.rules.items() if t in set(mod.genes)})
        results["claude"] = evaluate_spec(mod, spec, "claude")
        out["specs"]["claude"] = json.loads(spec.to_json())
        print(f"  claude done ({len(spec.edges)} edges)")
        flush()

    need_pool = ("oracle" in sources) or ("random" in sources)
    if need_pool:
        # spawn, not fork. Evaluating any structural source above initialises JAX in
        # this process, and forking a process with XLA's thread pool live deadlocks
        # the children on their first compile. Spawned workers import cleanly.
        ctx = mp.get_context("spawn")
        with ctx.Pool(args.workers, initializer=_init_worker,
                      initargs=(mod.genes, mod.perts, mod.observed, mod.de_mask)) as pool:

            if "oracle" in sources:
                print(f"  oracle: {args.seeds} seeds on {args.workers} workers")
                seeds_out, best_overall, best_score = [], None, -np.inf
                for s in range(args.seeds):
                    print(f"  seed {s}")
                    res, trace = run_oracle(pool, mod, seed=s)
                    spec = oracle_search.spec_from_edges(mod.genes, res["edges"])
                    ev = evaluate_spec(mod, spec, f"oracle_seed{s}")
                    final = compare.bootstrap_mean(ev.scores["de_overlap"])
                    res["final_de_overlap"] = final
                    seeds_out.append(res)
                    if final["mean"] > best_score:
                        best_score, best_overall = final["mean"], (res, spec, ev, trace)
                    print(f"    final LOO DE-overlap {final['mean']:.4f} "
                          f"[{final['lo']:.4f}, {final['hi']:.4f}]")
                    out["oracle_seeds"] = seeds_out
                    flush()

                res, spec, ev, trace = best_overall
                results["oracle"] = compare.SourceResult(
                    "oracle", ev.predictions, ev.scores, n_edges=len(spec.edges),
                    meta={"anneal_gain": round(res["anneal_gain"], 4)})
                out["specs"]["oracle"] = json.loads(spec.to_json())
                out["oracle_ceiling_spread"] = {
                    "by_seed": [r["final_de_overlap"]["mean"] for r in seeds_out],
                    "sd": float(np.std([r["final_de_overlap"]["mean"] for r in seeds_out])),
                    "min": float(np.min([r["final_de_overlap"]["mean"] for r in seeds_out])),
                    "max": float(np.max([r["final_de_overlap"]["mean"] for r in seeds_out])),
                }
                out["diagnostics"]["from_trace"] = {
                    "equivalence_width": identifiability.equivalence_width(
                        trace.train_loss, trace.holdout),
                    "sign_stability": identifiability.sign_stability(
                        trace.edge_sets, trace.train_loss),
                }
                flush()

                # nested honest variant: select on an inner split, score on an outer one
                n = len(mod.perts)
                rng = np.random.default_rng(12345)
                order = rng.permutation(n)
                n_outer = max(2, int(round(NESTED_OUTER_FRACTION * n)))
                outer, inner = sorted(order[:n_outer].tolist()), sorted(order[n_outer:].tolist())
                print(f"  nested oracle: inner {len(inner)}, outer {len(outer)}")
                nres, _ = run_oracle(pool, mod, seed=0, subset=inner)
                nspec = oracle_search.spec_from_edges(mod.genes, nres["edges"])
                nev = evaluate_spec(mod, nspec, "oracle_nested")
                outer_scores = nev.scores["de_overlap"][outer]
                out["oracle_nested"] = {
                    "inner_perts": [mod.perts[i] for i in inner],
                    "outer_perts": [mod.perts[i] for i in outer],
                    "n_edges": len(nspec.edges),
                    "outer_de_overlap": compare.bootstrap_mean(outer_scores),
                }
                out["specs"]["oracle_nested"] = json.loads(nspec.to_json())
                flush()

            if "random" in sources:
                n_edges = args.random_edges or (
                    results["oracle"].n_edges if "oracle" in results else 24)
                print(f"  random null: {args.n_random} structures at {n_edges} edges")
                rng = np.random.default_rng(999)
                specs = []
                for _ in range(args.n_random):
                    try:
                        specs.append(random_null.sample_spec(mod.genes, mod.perts, n_edges, rng))
                    except ValueError as e:
                        print(f"    stopped sampling: {e}")
                        break
                tasks = [(tuple(oracle_search.edges_of(s)), len(mod.perts),
                          compare.FIT_STARTS, compare.FIT_MAX_ITER, 0, None) for s in specs]
                t0 = time.time()
                scored = pool.map(_score_edge_set, tasks, chunksize=1)
                vals = [h for h, _ in scored]
                out["random_null"] = {
                    "n_edges": n_edges,
                    "n_structures": len(vals),
                    "seconds": round(time.time() - t0, 1),
                    **random_null.summarise_null(vals),
                    "values": [round(v, 4) if v == v else None for v in vals],
                }
                print(f"  random null done in {time.time() - t0:.0f}s: "
                      f"mean {out['random_null'].get('mean', float('nan')):.4f}")
                flush()

    # comparator table and the pairwise structure agreement
    if "linear" in results:
        out["table"] = compare.comparator_table(results, reference="linear")
        out["table_vs_mean"] = compare.comparator_table(results, reference="mean") \
            if "mean" in results else None
        if "random_null" in out:
            for row in out["table"]:
                row["random_null_percentile"] = round(random_null.percentile_of(
                    row["de_overlap_mean"], out["random_null"]["values"]), 1)

    specs_for_jaccard = {k: ModelSpec.model_validate(v) for k, v in out["specs"].items()}
    out["structure_agreement"] = {
        f"{a}|{b}": round(compare.jaccard_edges(specs_for_jaccard[a], specs_for_jaccard[b]), 4)
        for i, a in enumerate(specs_for_jaccard) for b in list(specs_for_jaccard)[i + 1:]
    }
    out["per_fold_scores"] = {k: [None if v != v else round(float(v), 4)
                                  for v in r.scores["de_overlap"]]
                              for k, r in results.items()}
    flush()

    print("\n  source                de_overlap [95% CI]        delta vs linear [95% CI]   adv")
    for row in out.get("table", []):
        print(f"  {row['source']:<20} {row['de_overlap_mean']:.4f} "
              f"[{row['de_overlap_lo']:.4f}, {row['de_overlap_hi']:.4f}]   "
              f"{row['delta_vs_linear']:+.4f} [{row['delta_lo']:+.4f}, {row['delta_hi']:+.4f}]  "
              f"{'YES' if row['advantage'] else 'no'}")
    print(f"\nwrote {outpath}")


if __name__ == "__main__":
    main()
