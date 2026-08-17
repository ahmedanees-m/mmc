"""Step 5: is proposed structure conditioned on the data, or recalled from the names?

PREREG_v4 section 6. Four arms, one harness, the same modules and seeds, changing
exactly one thing at a time:

    A1  real symbols      real data          the status quo
    A2  aliased symbols   real data          does the proposal survive without names
    A3  real symbols      permuted data      does it change when the data is destroyed
    A4  aliased symbols   permuted data      floor

The comparison that carries the claim is between-arm edge agreement measured against
the within-arm ceiling. J(A1, A1'), where A1' is A1 at another seed, is how much
agreement an unperturbed replicate produces. J(A1, A2) below that means the names are
doing work; J(A1, A3) at that ceiling means the data is not.

Before any arm runs, the assembled prompt for A2 and A4 is audited for surviving gene
symbols and the run aborts if any are found. An anonymisation arm that leaks is worse
than no arm, because it produces a confident null.

    python scripts/step5_anonymization.py --modules Cytokine_production,TCR_signalosome \\
        --seeds 3 --out results/step5.json

Requires MMC_ZHU_STORE and a model key. Every call is added to paper/LEDGER_api_spend.md.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np

from mmc.data import module_data, module_extract
from mmc.eval import compare
from mmc.loop import anonymize as anon
from mmc.loop.run import discover
from mmc.shared import llm

ARMS = ("A1", "A2", "A3", "A4")
ARM_SPEC = {"A1": (False, False), "A2": (True, False),
            "A3": (False, True), "A4": (True, True)}


def _observed_from_matrix(mod, genes, obs) -> dict:
    """The nested-dict shape the loop consumes, from an aligned matrix."""
    return {
        mod.perts[i]: {genes[j]: float(obs[i, j])
                       for j in range(len(genes)) if genes[j] != mod.perts[i]}
        for i in range(len(mod.perts))
    }


def build_arm(mod, module: str, arm: str, seed: int):
    """Return (genes, observed, alias_map, audit) for one arm."""
    aliased, permuted = ARM_SPEC[arm]
    obs = np.asarray(mod.observed, float)
    de = np.asarray(mod.de_mask, bool)
    if permuted:
        obs, de = anon.permute_perturbation_labels(obs, de, seed=seed)

    amap = anon.make_alias_map(mod.genes, seed=seed) if aliased else None
    context = f"Module {module} in CD4+ T cells."
    if aliased:
        genes = [amap.to_alias[g] for g in mod.genes]
        perts = [amap.to_alias[p] for p in mod.perts]
        context = ("A set of genes measured in a perturbation experiment. Gene "
                   "identities are withheld.")
    else:
        genes, perts = list(mod.genes), list(mod.perts)

    observed = {
        perts[i]: {genes[j]: float(obs[i, j]) for j in range(len(genes))
                   if genes[j] != perts[i]}
        for i in range(len(perts))
    }
    surface = anon.assemble_prompt_surface(genes, context, None, "")
    audit = anon.redaction_violations(surface, list(mod.genes)) if aliased else []
    return genes, observed, context, amap, audit, de


def run_arm(module: str, condition: str, mod, arm: str, seed: int, *,
            max_iters: int, n_starts: int, max_iter: int) -> dict:
    genes, observed, context, amap, audit, de = build_arm(mod, module, arm, seed)
    if audit:
        raise RuntimeError(
            f"{module}/{arm}/seed{seed}: gene symbols survived redaction: {audit}")

    t0 = time.time()
    result = discover(module, condition, context=context, observed=observed,
                      genes=genes, max_iters=max_iters, n_starts=n_starts,
                      max_iter=max_iter)
    spec = result.ensemble.best().spec
    if amap is not None:
        spec = anon.deanonymise_spec(spec, amap)

    # every arm is scored against the real, unpermuted data, so the arms differ only
    # in what the proposer saw and not in what its structure is judged on
    fn = compare.bind_structural(spec, mod)
    res = compare.evaluate_source(mod, arm, fn, n_edges=len(spec.edges))
    marg = compare.bootstrap_mean(res.scores["de_overlap"])
    return {
        "arm": arm, "seed": seed, "module": module,
        "n_edges": len(spec.edges),
        "edges": sorted((e.regulator, e.target, e.sign) for e in spec.edges),
        "de_overlap": marg,
        "per_fold": [None if v != v else round(float(v), 4)
                     for v in res.scores["de_overlap"]],
        "iterations": len(result.history),
        "seconds": round(time.time() - t0, 1),
        "redaction_clean": True,
    }


def jaccard(a: list, b: list) -> float:
    sa = {(r, t) for r, t, _ in a}
    sb = {(r, t) for r, t, _ in b}
    if not sa and not sb:
        return float("nan")
    return len(sa & sb) / len(sa | sb)


def agreement_matrix(runs: list[dict]) -> dict:
    """Between-arm and within-arm edge agreement, the statistic that carries C5."""
    by = {}
    for r in runs:
        by.setdefault((r["module"], r["arm"]), {})[r["seed"]] = r["edges"]

    within, between = {}, {}
    for (module, arm), seeds in by.items():
        vals = [jaccard(a, b) for s1, s2 in itertools.combinations(sorted(seeds), 2)
                for a, b in [(seeds[s1], seeds[s2])]]
        vals = [v for v in vals if v == v]
        if vals:
            within[f"{module}/{arm}"] = {"mean": float(np.mean(vals)), "n_pairs": len(vals)}

    modules = {m for m, _ in by}
    for module in modules:
        for arm in ARMS[1:]:
            pairs = []
            for seed, edges in by.get((module, "A1"), {}).items():
                other = by.get((module, arm), {}).get(seed)
                if other is not None:
                    pairs.append(jaccard(edges, other))
            pairs = [v for v in pairs if v == v]
            if pairs:
                between[f"{module}/A1_vs_{arm}"] = {"mean": float(np.mean(pairs)),
                                                    "n_seeds": len(pairs)}
    return {"within_arm_seed_ceiling": within, "between_arm": between}


def interpret(agree: dict) -> list[str]:
    """State the reading the pre-registration attaches to each pattern."""
    notes = []
    for key, val in sorted(agree["between_arm"].items()):
        module = key.split("/")[0]
        ceiling = agree["within_arm_seed_ceiling"].get(f"{module}/A1", {}).get("mean")
        if ceiling is None:
            continue
        ratio = val["mean"] / ceiling if ceiling else float("nan")
        notes.append(f"{key}: {val['mean']:.3f} against a within-A1 ceiling of "
                     f"{ceiling:.3f} (ratio {ratio:.2f})")
    return notes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--modules", required=True)
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--max-iters", type=int, default=2)
    ap.add_argument("--n-starts", type=int, default=4)
    ap.add_argument("--max-iter", type=int, default=200)
    ap.add_argument("--module-def", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(d["module"], d["regulators"], d["targets"])

    modules = [m.strip() for m in args.modules.split(",") if m.strip()]
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    out = {"condition": args.condition, "arms": arms, "seeds": args.seeds,
           "model": llm.model(), "runs": [], "errors": []}
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    def flush():
        with open(outpath, "w") as f:
            json.dump(out, f, indent=2, default=float)

    for module in modules:
        mod = module_data.build_module_data(module, args.condition)
        print(f"=== {module}: {len(mod.genes)} genes, {len(mod.perts)} perturbations ===")
        for seed in range(args.seeds):
            for arm in arms:
                try:
                    r = run_arm(module, args.condition, mod, arm, seed,
                                max_iters=args.max_iters, n_starts=args.n_starts,
                                max_iter=args.max_iter)
                except Exception as e:                          # noqa: BLE001
                    print(f"  {arm} seed {seed}: FAILED {e}")
                    out["errors"].append({"module": module, "arm": arm, "seed": seed,
                                          "error": str(e)})
                    flush()
                    continue
                out["runs"].append(r)
                print(f"  {arm} seed {seed}: {r['n_edges']} edges, held-out "
                      f"{r['de_overlap']['mean']:.4f} "
                      f"[{r['de_overlap']['lo']:.4f}, {r['de_overlap']['hi']:.4f}], "
                      f"{r['seconds']:.0f}s")
                flush()

    if out["runs"]:
        out["agreement"] = agreement_matrix(out["runs"])
        out["reading"] = interpret(out["agreement"])
        flush()
        print("\nedge agreement:")
        for line in out["reading"]:
            print("  " + line)
    print(f"\nwrote {outpath}")


if __name__ == "__main__":
    main()
