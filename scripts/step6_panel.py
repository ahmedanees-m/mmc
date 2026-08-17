"""Step 6: is grounded-but-non-predictive a property of one model or of the class?

PREREG_v4 section 7. The reliability claim currently rests on one model family. This
runs the identical loop, prompt, schema validation and retry path across a panel of
model families and reports per-model and pooled statistics.

The identifiers are locked in amendment A6 against verified availability, not against
the catalogue: of 41 plausible candidates probed, 18 served a completion and 23 did
not. No DeepSeek model is served on this key, so the reasoning slot substitutes an
NVIDIA Nemotron model, which section 7 permits provided the substitution is recorded.

One prompt for every model. Section 7 forbids per-model tuning, and there is
deliberately nowhere in this script to put it.

The metric that carries the claim is the held-out validated rate: the fraction of
proposals whose model beats the linear baseline on held-out DE-overlap by the section
1.2 rule. Under Branch A that rate is expected to be zero for every model, and the
value of the panel is the tightened pooled interval and whether the models order
themselves by how much of the oracle ceiling they attain.

    python scripts/step6_panel.py --modules Cytokine_production --seeds 2 \\
        --ceiling results/step1_Cytokine_production.json --out results/step6.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

import numpy as np

from mmc.data import module_data, module_extract
from mmc.eval import compare, random_null
from mmc.loop import providers
from mmc.loop.run import discover

# Amendment A6. Anthropic tiers first, then the NGC panel.
PANEL = [
    ("claude-opus-4-8", "anthropic"),
    ("claude-sonnet-5", "anthropic"),
    ("openai/gpt-oss-120b", "ngc"),
    ("meta/llama-3.1-70b-instruct", "ngc"),
    ("nvidia/nemotron-3-ultra-550b-a55b", "ngc"),
    ("z-ai/glm-5.2", "ngc"),
]


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval, the form used for the corpus rate throughout."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def run_one(module: str, condition: str, mod, model: str, provider: str, seed: int,
            linear_scores: np.ndarray, *, max_iters: int, n_starts: int,
            max_iter: int) -> dict:
    os.environ["MMC_MODEL"] = model
    os.environ["MMC_PROVIDER"] = provider
    before = providers.STATS.as_dict()

    t0 = time.time()
    result = discover(module, condition, max_iters=max_iters, n_starts=n_starts,
                      max_iter=max_iter)
    spec = result.ensemble.best().spec

    fn = compare.bind_structural(spec, mod)
    res = compare.evaluate_source(mod, model, fn, n_edges=len(spec.edges))
    scores = res.scores["de_overlap"]
    delta = compare.paired_delta(scores, linear_scores)
    after = providers.STATS.as_dict()

    return {
        "model": model, "provider": provider, "module": module, "seed": seed,
        "n_edges": len(spec.edges),
        "edges": sorted((e.regulator, e.target, e.sign) for e in spec.edges),
        "de_overlap": compare.bootstrap_mean(scores),
        "delta_vs_linear": delta,
        "validated": bool(delta["advantage"]),
        "iterations": len(result.history),
        "calls": after["calls"] - before["calls"],
        "transport_retries": after["transport_retries"] - before["transport_retries"],
        "seconds": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--modules", required=True)
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--models", default="")
    ap.add_argument("--ceiling", default="", help="a Step 1 result JSON, for the "
                                                 "fraction-of-ceiling column")
    ap.add_argument("--max-iters", type=int, default=2)
    ap.add_argument("--n-starts", type=int, default=4)
    ap.add_argument("--max-iter", type=int, default=200)
    ap.add_argument("--random-control", type=int, default=50)
    ap.add_argument("--module-def", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(d["module"], d["regulators"], d["targets"])

    panel = PANEL
    if args.models:
        wanted = {m.strip() for m in args.models.split(",")}
        panel = [(m, p) for m, p in PANEL if m in wanted]

    ceilings = {}
    if args.ceiling:
        with open(args.ceiling) as f:
            c = json.load(f)
        for row in c.get("table", []):
            if row["source"] == "oracle":
                ceilings[c["module"]] = row["de_overlap_mean"]

    modules = [m.strip() for m in args.modules.split(",") if m.strip()]
    out = {"condition": args.condition, "panel": [m for m, _ in panel],
           "seeds": args.seeds, "runs": [], "errors": [], "ceilings": ceilings}
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    def flush():
        with open(outpath, "w") as f:
            json.dump(out, f, indent=2, default=float)

    for module in modules:
        mod = module_data.build_module_data(module, args.condition)
        linear = compare.evaluate_source(mod, "linear", compare.bind_linear(mod))
        lin = linear.scores["de_overlap"]
        out.setdefault("linear", {})[module] = compare.bootstrap_mean(lin)

        # the reasoning-versus-search control: random structures at the panel's own
        # edge counts, so a model is not credited for beating nothing
        rng = np.random.default_rng(0)
        ctrl = []
        for _ in range(args.random_control):
            spec = random_null.sample_spec(mod.genes, mod.perts, 20, rng)
            fn = compare.bind_structural(spec, mod)
            ctrl.append(compare.bootstrap_mean(
                compare.score_predictions(mod, compare.fold_predictions(mod, fn))
                ["de_overlap"])["mean"])
        out.setdefault("random_control", {})[module] = random_null.summarise_null(ctrl)
        print(f"=== {module}: linear {out['linear'][module]['mean']:.4f}, "
              f"random control {out['random_control'][module].get('mean', float('nan')):.4f} ===")
        flush()

        for model, provider in panel:
            for seed in range(args.seeds):
                try:
                    r = run_one(module, args.condition, mod, model, provider, seed,
                                lin, max_iters=args.max_iters, n_starts=args.n_starts,
                                max_iter=args.max_iter)
                except Exception as e:                          # noqa: BLE001
                    print(f"  {model} seed {seed}: FAILED {type(e).__name__}: {e}")
                    out["errors"].append({"model": model, "module": module,
                                          "seed": seed, "error": str(e)[:300]})
                    flush()
                    continue
                out["runs"].append(r)
                d = r["delta_vs_linear"]
                print(f"  {model:<38} seed {seed}: {r['n_edges']:>3} edges  "
                      f"held-out {r['de_overlap']['mean']:.4f}  "
                      f"delta {d['delta']:+.4f} [{d['lo']:+.4f}, {d['hi']:+.4f}]  "
                      f"{'VALIDATED' if r['validated'] else 'not validated'}  "
                      f"{r['seconds']:.0f}s")
                flush()

    # per-model and pooled statistics
    per_model = {}
    for model, _ in panel:
        runs = [r for r in out["runs"] if r["model"] == model]
        if not runs:
            per_model[model] = {"n_runs": 0, "attempted": True, "served": False}
            continue
        k = sum(1 for r in runs if r["validated"])
        lo, hi = wilson(k, len(runs))
        frac = []
        for r in runs:
            ceil = ceilings.get(r["module"])
            if ceil:
                frac.append(r["de_overlap"]["mean"] / ceil)
        per_model[model] = {
            "n_runs": len(runs), "served": True,
            "validated": k, "validated_rate": k / len(runs),
            "validated_wilson95": [lo, hi],
            "mean_de_overlap": float(np.mean([r["de_overlap"]["mean"] for r in runs])),
            "mean_edges": float(np.mean([r["n_edges"] for r in runs])),
            "fraction_of_ceiling": float(np.mean(frac)) if frac else None,
            "calls": sum(r["calls"] for r in runs),
            "transport_retries": sum(r["transport_retries"] for r in runs),
        }
    out["per_model"] = per_model

    total = len(out["runs"])
    validated = sum(1 for r in out["runs"] if r["validated"])
    lo, hi = wilson(validated, total)
    out["pooled"] = {"n_runs": total, "validated": validated,
                     "validated_rate": (validated / total) if total else float("nan"),
                     "wilson95": [lo, hi],
                     "call_stats": providers.STATS.as_dict()}
    flush()

    print(f"\n{'model':<40}{'runs':>5}{'valid':>7}{'held-out':>10}{'frac ceiling':>14}"
          f"{'retries':>9}")
    for model, e in per_model.items():
        if not e.get("served"):
            print(f"  {model:<38}{'-':>5}  did not serve a completion")
            continue
        fc = e["fraction_of_ceiling"]
        print(f"  {model:<38}{e['n_runs']:>5}{e['validated']:>7}"
              f"{e['mean_de_overlap']:>10.4f}"
              f"{(f'{fc:.2f}' if fc else 'n/a'):>14}{e['transport_retries']:>9}")
    p = out["pooled"]
    print(f"\npooled validated rate {p['validated']}/{p['n_runs']} "
          f"= {p['validated_rate']:.3f}, Wilson 95% "
          f"[{p['wilson95'][0]:.3f}, {p['wilson95'][1]:.3f}]")
    print(f"wrote {outpath}")


if __name__ == "__main__":
    main()
