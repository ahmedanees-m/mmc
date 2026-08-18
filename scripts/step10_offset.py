"""Step 10, promoted: does structure earn its keep once the shared response is handed to it?

Section 2.2 found that a sparse structural model cannot represent the stereotyped bulk
response at all, because its prediction is the difference between clamped and unclamped
fixed points. The mean baseline is that bulk response and nothing else, and it reaches a
large fraction of what the linear map achieves. So the structural arm forfeits the
dominant component of the signal before any search begins.

`structural + offset` hands that component over: the same structure and the same solve,
plus a per-gene offset fitted on the training folds of each split. If it closes the gap
to linear, the ceiling was the forfeit and the grammar is not the limitation. If it does
not, the shortfall lies elsewhere and the architectural explanation is wrong.

Structures are re-used from the recorded comparator results rather than searched again,
so this adds no search and the structures are exactly the ones already reported.

    python scripts/step10_offset.py --result results/step1_Cytokine_production_full.json \
        --sources claude,oracle --out results/step10_cytokine.json
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
    ap.add_argument("--module-def", default="")
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
        # modules registered dynamically for the original run are not in the static
        # registry, but every recorded structure carries the gene list it was built over,
        # so the module is rebuilt from the result itself rather than re-derived
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
    # the baselines this has to be read against, recomputed here so every number in the
    # table comes from one run rather than being carried across from another
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
