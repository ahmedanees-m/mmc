"""Count the modules where per-fold structure selection exceeds the linear map.

The per-fold arm is scored on its outer split, while the linear arm is recorded under
leave-one-out over every fold. Comparing the two directly pairs different fold sets, so
three pairings are reported and the difference between them is what the script shows:

  outer against all      the per-fold arm on its outer split against the linear arm on
                         every fold, which pairs different fold sets
  all against all        the selected structure re-scored on every fold, matching the
                         linear arm's folds. Matched, but the selection saw the inner
                         perturbations, so scoring it on inner folds is partly leaky.
  outer against outer    both arms restricted to the outer split. Matched, and on the
                         folds the selection never saw.

No refitting is involved. The linear arm's per-fold scores are already recorded under
leave-one-out, so each fold's value is a prediction made without that perturbation, and
taking the subset of folds belonging to the outer split is a selection over recorded
numbers.

The outer splits are small, so the third pairing carries little power per module. The
fold-level pooled estimate over the same splits is in pool_nested_folds.py.

    python scripts/run_nested_search.py --dir results/comparator/twenty_seeds \
        --offset-dir results/comparator/offset --out results/nested/nested_matched.json

Reads result JSON only; no model is refit and no data store is needed.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/comparator/twenty_seeds")
    ap.add_argument("--offset-dir", default="results/comparator/offset")
    ap.add_argument("--out", default="results/nested/nested_matched.json")
    args = ap.parse_args()

    rows = []
    for p in sorted(glob.glob(f"{args.dir}/*.json")):
        d = json.load(open(p))
        rn = d["random_null"]
        table = {r["source"]: r for r in d["table"]}
        lin_all = table["linear"]["de_overlap_mean"]
        if lin_all <= rn["p95"]:
            continue

        nested = d["oracle_nested"]
        perts = d["module_summary"]["perts"]
        idx = {g: i for i, g in enumerate(perts)}
        outer = [idx[g] for g in nested["outer_perts"] if g in idx]

        lin_folds = d["per_fold_scores"]["linear"]
        lin_outer = [lin_folds[i] for i in outer
                     if i < len(lin_folds) and lin_folds[i] is not None]

        # the selected structure re-scored on the full fold set, from the offset run
        off = Path(args.offset_dir) / f"real_{d['module']}.json"
        nested_all = None
        if off.exists():
            for r in json.load(open(off))["table"]:
                if r["source"] == "oracle_nested":
                    nested_all = r["de_overlap_mean"]

        rows.append({
            "module": d["module"],
            "n_outer_scoreable": len(lin_outer),
            "n_outer_perts": len(outer),
            "linear_all_folds": lin_all,
            "linear_outer_folds": float(np.mean(lin_outer)) if lin_outer else None,
            "nested_outer_folds": nested["outer_de_overlap"]["mean"],
            "nested_all_folds": nested_all,
        })

    def count(a, b):
        return sum(1 for r in rows
                   if r[a] is not None and r[b] is not None and r[a] > r[b])

    print(f"{'module':<24}{'outer n':>8}{'linear all':>12}{'linear outer':>14}"
          f"{'nested outer':>14}{'nested all':>12}")
    for r in rows:
        lo = f"{r['linear_outer_folds']:.4f}" if r["linear_outer_folds"] is not None else "n/a"
        na = f"{r['nested_all_folds']:.4f}" if r["nested_all_folds"] is not None else "n/a"
        print(f"{r['module']:<24}{r['n_outer_scoreable']:>8}{r['linear_all_folds']:>12.4f}"
              f"{lo:>14}{r['nested_outer_folds']:>14.4f}{na:>12}")

    n = len(rows)
    mismatched = count("nested_outer_folds", "linear_all_folds")
    all_folds = count("nested_all_folds", "linear_all_folds")
    matched = count("nested_outer_folds", "linear_outer_folds")

    print(f"\nper-fold arm exceeds linear, {n} sound modules:")
    print(f"  outer folds against all folds   {mismatched} of {n}   different fold sets")
    print(f"  all folds against all folds     {all_folds} of {n}   matched, partly leaky")
    print(f"  outer folds against outer folds {matched} of {n}   matched and out of fold")

    scoreable = [r["n_outer_scoreable"] for r in rows]
    print(f"\nscoreable outer folds per module: {sorted(scoreable)}, median "
          f"{np.median(scoreable):g}. The third pairing rests on those, so it is like for "
          f"like\nbut carries little power per module.")

    out = {"n_modules": n, "counts": {"nested_outer_vs_linear_all": mismatched,
                                      "nested_all_vs_linear_all": all_folds,
                                      "nested_outer_vs_linear_outer": matched},
           "rows": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
