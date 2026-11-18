"""Re-score recorded structures with a fitted per-gene offset added to the prediction.

A prediction in this model class is the difference between a clamped and an unclamped fixed
point, so it cannot express a response shared across perturbations. The mean baseline is that
shared response and nothing else, and on the primary module it reaches 78% of what the linear
map achieves. The structural arm therefore forfeits the dominant component of the signal before
any structure is chosen.

This script hands that component back. The same structure and the same solve are used, plus a
per-gene offset fitted on the training rows of each fold, which is the one quantity the grammar
cannot emit and nothing else. If the offset closes the gap to the linear map, the forfeit
accounts for the gap; if it does not, the shortfall lies elsewhere.

Structures are re-used from a recorded comparator result rather than searched again, so this
repeats no search and the structures are exactly the ones already reported.

**Every arm here is scored on the module's full fold set.** For `oracle` that is the fold set it
was selected on, and for `oracle_nested` it includes the inner folds its selection used, so these
numbers are not out-of-fold estimates and are not interchangeable with the outer-fold numbers in
`results/nested/`. What the comparison isolates is the size of the offset's contribution with the
structure held fixed, which is the same for every arm because the fold set is.

    python scripts/run_offset_arm.py \
        --result results/comparator/twenty_seeds/Cytokine_production.json \
        --sources oracle,oracle_nested \
        --out results/comparator/offset/real_Cytokine_production.json

Reads per-module response matrices, either from the exported matrices named by
MMC_MODULE_MATRICES or from the store named by MMC_ZHU_STORE. Structures come from the result
file, so this script makes no API calls and is deterministic.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from mmc.data import module_data, module_extract
from mmc.eval import compare, model_classes
from mmc.grammar.model_spec import ModelSpec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", required=True, help="comparator result carrying the specs")
    ap.add_argument("--sources", default="oracle",
                    help="comma-separated spec names to re-evaluate with an offset")
    ap.add_argument("--module-def", default="",
                    help="JSON registering a module that is not in the static registry")
    ap.add_argument("--condition", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    prior = json.load(open(args.result))
    module = prior["module"]
    condition = args.condition or prior.get("condition", "Stim8hr")

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(module, d["regulators"], d["targets"])
    elif module not in module_extract.MODULES:
        # A module assembled from the atlas is not in the static registry, but every
        # recorded structure carries the gene list it was built over, so the module is
        # rebuilt from the result file itself rather than re-derived.
        genes = next((s["genes"] for s in (prior.get("specs") or {}).values()
                      if isinstance(s, dict) and s.get("genes")), None)
        if genes is None:
            raise SystemExit(f"{module} is not registered and {args.result} records no "
                             f"gene list to rebuild it from")
        module_extract.register_module(module, genes, genes)
        print(f"registered {module} from the recorded structure: {len(genes)} genes")

    mod = module_data.build_module_data(module, condition)
    summary = module_data.module_summary(mod)
    print(f"{module} / {condition}: {summary['n_genes']} genes, {summary['n_perts']} "
          f"perturbations, {summary['n_folds_with_de']} scoreable folds")

    results = {}
    # The baselines this has to be read against are recomputed here, so every number in the
    # table comes from one run rather than being carried across from another.
    results["linear"] = compare.evaluate_source(mod, "linear", compare.bind_linear(mod))
    results["mean"] = compare.evaluate_source(mod, "mean", compare.bind_mean(mod))
    results["zero"] = compare.evaluate_source(mod, "zero", compare.bind_zero(mod))

    specs = prior.get("specs") or {}
    for name in [s.strip() for s in args.sources.split(",") if s.strip()]:
        if name not in specs:
            print(f"  {name}: no structure recorded in {args.result}, skipped")
            continue
        spec = ModelSpec.from_json(json.dumps(specs[name]))
        n_edges = len(spec.edges)
        results[name] = compare.evaluate_source(
            mod, name, compare.bind_structural(spec, mod), n_edges=n_edges)
        results[f"{name}+offset"] = compare.evaluate_source(
            mod, f"{name}+offset",
            model_classes.bind_structural_with_offset(spec, mod), n_edges=n_edges)
        print(f"  {name}: {n_edges} edges, evaluated with and without the offset")

    table = compare.comparator_table(results, reference="linear")
    vs_mean = compare.comparator_table(results, reference="mean")
    out = {"module": module, "condition": condition, "module_summary": summary,
           "table": table, "table_vs_mean": vs_mean,
           "per_fold_scores": {k: list(v.scores["de_overlap"]) for k, v in results.items()}}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=2, default=float)

    print(f"\n{'source':<22}{'DE-overlap':>12}{'vs linear':>24}{'adv':>6}")
    for r in sorted(table, key=lambda r: -r["de_overlap_mean"]):
        print(f"  {r['source']:<20}{r['de_overlap_mean']:>12.4f}"
              f"   {r['delta_vs_linear']:+.4f} [{r['delta_lo']:+.4f}, {r['delta_hi']:+.4f}]"
              f"{'  YES' if r['advantage'] else '   no':>6}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
