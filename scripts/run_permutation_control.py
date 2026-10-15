"""Amendment A22: does the oracle ceiling survive destroying the signal?

The leaky oracle selects its structure on the folds it is then scored on, over a search that
evaluates a median of 13,271 structures at the full objective. Its advantage over a linear
map therefore bundles two things that no previous analysis separated: real structure in the
data, and the maximum of a large selection made on the scored folds.

The permuted arm runs the identical search, at the same budget, over the same structure
space, on data whose perturbation-to-response pairing has been destroyed. Everything that
could produce an advantage by selection survives; only the signal is removed. Whatever the
permuted search reaches is what selection alone buys.

Two questions are answered separately, and they are not the same question:

  does the real search beat the permuted one?   whether any signal is being used at all
  does the PERMUTED search beat the linear map? whether "the ceiling exceeds linear" is
                                                evidence of recoverable structure

The second is the one claim 1 rests on. If a search with no signal also clears the linear
map, then clearing the linear map is a property of the search rather than of the data.

    python scripts/step_permutation_verdict.py --out results/a22_verdict.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", default="results/a21")
    ap.add_argument("--permuted", default="results/a22")
    ap.add_argument("--out", default="results/a22_verdict.json")
    args = ap.parse_args()

    rows = []
    for pf in sorted(glob.glob(f"{args.permuted}/step1_*.json")):
        perm = json.load(open(pf))
        module = perm["module"]
        rf = os.path.join(args.real, os.path.basename(pf))
        if not os.path.exists(rf):
            continue
        real = json.load(open(rf))
        lin = {r["source"]: r for r in real["table"]}["linear"]["de_overlap_mean"]
        if lin <= real["random_null"]["p95"]:
            continue

        real_max = max(s["final_de_overlap"]["mean"] for s in real["oracle_seeds"])
        perm_seeds = [s["final_de_overlap"]["mean"] for s in perm["oracle_seeds"]]
        perm_max = max(perm_seeds)
        perm_lin = {r["source"]: r for r in perm["table"]}["linear"]["de_overlap_mean"]

        rows.append({
            "module": module,
            "linear": lin,
            "real_oracle_max": real_max,
            "permuted_oracle_max": perm_max,
            "permuted_seeds": [round(v, 4) for v in perm_seeds],
            "permuted_linear": perm_lin,
            "ratio_permuted_to_real": perm_max / real_max if real_max else float("nan"),
            "real_beats_permuted": bool(real_max > perm_max),
            "permuted_beats_linear": bool(perm_max > lin),
            "real_minus_permuted": real_max - perm_max,
        })

    n = len(rows)
    print(f"{'module':<24}{'linear':>9}{'real max':>10}{'permuted':>10}{'ratio':>8}"
          f"{'real-perm':>11}  perm>linear")
    for r in sorted(rows, key=lambda x: -x["ratio_permuted_to_real"]):
        print(f"{r['module']:<24}{r['linear']:>9.4f}{r['real_oracle_max']:>10.4f}"
              f"{r['permuted_oracle_max']:>10.4f}{r['ratio_permuted_to_real']:>8.2f}"
              f"{r['real_minus_permuted']:>+11.4f}  {r['permuted_beats_linear']}")

    beats = sum(r["real_beats_permuted"] for r in rows)
    perm_clears = sum(r["permuted_beats_linear"] for r in rows)
    ratios = [r["ratio_permuted_to_real"] for r in rows]

    print(f"\nmodules: {n}")
    print(f"real search exceeds the permuted search:  {beats} of {n}")
    print(f"PERMUTED search exceeds the linear map:   {perm_clears} of {n}")
    print(f"permuted as a fraction of real: median {np.median(ratios):.2f}, "
          f"range {min(ratios):.2f} to {max(ratios):.2f}")
    print("\nThe permuted arm has no perturbation-to-response signal. Whatever it clears, it "
          "clears\nby selecting the best of a large search on the folds it is scored on. A "
          "count above zero\non the second line means the oracle exceeding the linear map is "
          "not by itself evidence\nthat recoverable structure is present.")
    n_seeds = len(rows[0]["permuted_seeds"]) if rows else 0
    if n_seeds and n_seeds < 20:
        print(f"\nThe permuted arm ran {n_seeds} seeds against the real arm's 20, so its "
              "maximum is drawn from a\nsmaller pool and is understated relative to a matched "
              "comparison.")
    else:
        print("\nBoth arms ran 20 seeds, so the maxima are drawn from equally sized pools and "
              "the\ncomparison is matched.")

    # The ratio of raw scores understates the point. What the paper claims is that the
    # signal-free search recovers most of the real search's ADVANTAGE OVER THE LINEAR MAP,
    # which is the quantity a reader cares about, so it is computed and stored here rather
    # than derived by hand at writing time.
    adv_real = [r["real_oracle_max"] - r["linear"] for r in rows]
    adv_perm = [r["permuted_oracle_max"] - r["linear"] for r in rows]
    frac = [b / a for a, b in zip(adv_real, adv_perm) if a > 0]
    pooled = sum(adv_perm) / sum(adv_real) if sum(adv_real) else float("nan")
    print(f"\nadvantage over the linear map: median {np.median(adv_real):+.4f} real, "
          f"{np.median(adv_perm):+.4f} signal-free")
    print(f"signal-free advantage as a fraction of real: median {np.median(frac):.2f} over the "
          f"{len(frac)} modules where the real advantage is positive, pooled {pooled:.2f}")

    out = {"n_modules": n, "n_real_beats_permuted": beats,
           "n_permuted_beats_linear": perm_clears,
           "median_ratio": float(np.median(ratios)),
           "median_advantage_real": float(np.median(adv_real)),
           "median_advantage_permuted": float(np.median(adv_perm)),
           "advantage_fraction_median": float(np.median(frac)),
           "advantage_fraction_pooled": float(pooled),
           "rows": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
