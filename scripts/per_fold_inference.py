"""Re-infer the two algorithmic structure sources inside each fold.

`scp_infer_adapter.algorithmic_edges` is normally called once on the whole module, so the
edge set is chosen using every perturbation including the one each fold holds out. The
parameters are then refit per fold, but the topology is not out of fold: `mean_difference`
ranks edges by the observed effect size, so for held-out perturbation i the edges out of
gene i were selected because of the very response the model is scored on predicting.

That matters because `mean_difference` is the one named structure source that exceeds the
linear map anywhere in this comparison. This script rebuilds both sources fold by fold from
the training rows only and reports the same comparison, with the whole-module form beside it.

    python scripts/per_fold_inference.py --module TCR_signalosome \
        --out results/comparator/per_fold_inference/TCR_signalosome.json

Reads per-module response matrices.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mmc.compare import scp_infer_adapter
from mmc.data import module_data
from mmc.eval import compare, oracle_search
from mmc.eval.holdout import ModuleData

MAX_EDGES = 30


def subset(mod: ModuleData, keep: list[int]) -> ModuleData:
    """A module holding only the given perturbation rows, for structure inference."""
    return ModuleData(genes=list(mod.genes), perts=[mod.perts[i] for i in keep],
                      observed=mod.observed[keep], de_mask=mod.de_mask[keep], spec=mod.spec)


def per_fold_source(mod: ModuleData, method: str):
    """Infer the structure on the training rows, fit on them, predict the held-out one."""
    from mmc.compile import structural
    from mmc.fit import fit_structural

    genes = list(mod.genes)
    edge_counts: list[int] = []

    def fn(train_idx, held_idx):
        train_mod = subset(mod, list(train_idx))
        edges = scp_infer_adapter.algorithmic_edges(train_mod, method, max_edges=MAX_EDGES)
        edge_counts.append(len(edges))
        if not edges:
            return np.zeros(len(genes), float)
        spec = oracle_search.spec_from_edges(genes, edges)
        observed = {
            mod.perts[i]: {genes[j]: float(mod.observed[i, j])
                           for j in range(len(genes)) if genes[j] != mod.perts[i]}
            for i in train_idx
        }
        fits = fit_structural.multi_fit(spec, observed, n_starts=compare.FIT_STARTS,
                                        max_iter=compare.FIT_MAX_ITER)
        return np.asarray(structural.knockdown(spec, fits[0]["params"], mod.perts[held_idx]))

    return fn, edge_counts


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--module-def", default="",
                    help="JSON registering a module, in the form run_comparator takes")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # Modules outside the static registry are passed the same way the comparator takes
    # them, so this runs on exactly the module the comparator ran.
    if args.module_def:
        from mmc.data import module_extract
        spec = json.load(open(args.module_def))
        module_extract.register_module(args.module, spec["regulators"], spec["targets"])

    mod = module_data.build_module_data(args.module, args.condition)
    print(f"=== per-fold inference on {args.module}/{args.condition}: "
          f"{len(mod.perts)} perturbations, {len(mod.genes)} genes ===", flush=True)

    linear = compare.evaluate_source(mod, "linear", compare.bind_linear(mod))
    out = {"module": args.module, "condition": args.condition,
           "linear": compare.bootstrap_mean(linear.scores["de_overlap"]),
           "per_fold_scores": {"linear": [None if v != v else round(float(v), 4)
                                          for v in linear.scores["de_overlap"]]},
           "sources": {}}

    for method in scp_infer_adapter.METHODS:
        # the whole-module form, exactly as the comparator runs it, so the two sit side by side
        leaky_edges = scp_infer_adapter.algorithmic_edges(mod, method, max_edges=MAX_EDGES)
        leaky = None
        if leaky_edges:
            spec = oracle_search.spec_from_edges(mod.genes, leaky_edges)
            res = compare.evaluate_source(mod, method, compare.bind_structural(spec, mod),
                                          n_edges=len(spec.edges))
            leaky = res.scores["de_overlap"]

        fn, edge_counts = per_fold_source(mod, method)
        scored = compare.evaluate_source(mod, f"{method}_per_fold", fn)
        h = scored.scores["de_overlap"]

        entry = {
            "n_edges_full_module": len(leaky_edges),
            "n_edges_per_fold": {"min": min(edge_counts), "max": max(edge_counts),
                                 "median": float(np.median(edge_counts))},
            "per_fold": compare.bootstrap_mean(h),
            "per_fold_vs_linear": compare.paired_delta(h, linear.scores["de_overlap"]),
            "per_fold_scores": [None if v != v else round(float(v), 4) for v in h],
        }
        if leaky is not None:
            entry["leaky"] = compare.bootstrap_mean(leaky)
            entry["leaky_vs_linear"] = compare.paired_delta(leaky,
                                                            linear.scores["de_overlap"])
            entry["leaky_minus_per_fold"] = compare.paired_delta(leaky, h)
        out["sources"][method] = entry

        lk = entry.get("leaky", {}).get("mean")
        print(f"  {method}: whole module {lk if lk is None else round(lk, 4)}  "
              f"per fold {round(entry['per_fold']['mean'], 4)}  "
              f"linear {round(out['linear']['mean'], 4)}", flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
