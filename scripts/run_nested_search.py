"""Resolve the nested-versus-linear count by scoring both arms on the same folds.

The record gives two different counts for the same comparison. Section 2.52 reports the
nested out-of-fold arm exceeding the linear map on 5 of 12 modules; the table gives 7 of
12. Neither is a mistake in arithmetic. They are two different comparisons:

  5 of 12   the nested arm scored on its outer split only, 1 to 8 folds with a median of 4,
            against the linear arm scored on the full leave-one-out folds, 7 to 19 of them.
            Different data, so the comparison is not like for like.
  7 of 12   the nested structure re-scored on the full leave-one-out folds, matching the
            linear arm's fold set. Matched, but the nested structure's search saw the inner
            perturbations, so scoring it on inner folds is partly leaky.

Neither is both out of fold and matched. This script produces the comparison that is: the linear
map restricted to each module's outer split, the same perturbations the nested arm was scored
on and the ones its search never saw.

No refitting is involved. The linear arm's per-fold scores are already recorded under
leave-one-out, so each fold's value is a prediction made without that perturbation; taking the
subset of folds belonging to the outer split is a selection over recorded numbers.

    python scripts/step_nested_matched.py --dir results/a21 --out results/nested_matched.json
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/a21")
    ap.add_argument("--offset-dir", default="results/a24")
    ap.add_argument("--out", default="results/nested_matched.json")
    args = ap.parse_args()

    rows = []
    for p in sorted(glob.glob(f"{args.dir}/step1_*.json")):
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

        # the nested structure re-scored on the full fold set, from the offset run
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

    print(f"\nnested exceeds linear, {n} sound modules:")
    print(f"  nested on outer folds vs linear on all folds   {mismatched} of {n}   "
          f"mismatched fold sets, the 2.52 figure")
    print(f"  nested on all folds   vs linear on all folds   {all_folds} of {n}   "
          f"matched folds, but partly leaky, the 2.53 figure")
    print(f"  nested on outer folds vs linear on outer folds {matched} of {n}   "
          f"**matched and honest**")

    scoreable = [r["n_outer_scoreable"] for r in rows]
    # Not int(): the median is 3.5, and truncating it to 3 is where the "median of three" in
    # the earlier record and draft came from.
    print(f"\nscoreable outer folds per module: {sorted(scoreable)}, median "
          f"{np.median(scoreable):g}. The matched comparison rests on those, so it is "
          f"honest\nand like for like but carries little power per module.")

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
