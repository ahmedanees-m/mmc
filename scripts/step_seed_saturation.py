"""Has the oracle search converged, and is its maximum over seeds a ceiling?

The fourth review makes the point that the seed maximum and the seed median answer different
questions. The maximum over k seeds is a downward-biased estimator of the true supremum, with
the bias shrinking as k grows, and it is the statistic the identifiability claim needs. The
median estimates the typical search outcome, which is a capability claim. Reporting both
leaves the reader to choose, and the choice is measurable rather than editorial.

The measurement: for each module and each k from 1 to the number of seeds, compute the mean
over all seed subsets of size k of the maximum within that subset. If the curve has flattened
by the largest k available, the maximum is a defensible estimate of the supremum. If it is
still climbing, the search has not converged and no ceiling has been estimated, in which case
the word has to come out of the paper and every claim rescopes to the structures actually
searched.

Uses the per-seed per-fold scores retained under amendment A20. No new compute.

    python scripts/step_seed_saturation.py --dir results/a20 --out results/seed_saturation.json
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/a20")
    ap.add_argument("--out", default="results/seed_saturation.json")
    args = ap.parse_args()

    rows = []
    for p in sorted(Path(args.dir).glob("step1_*.json")):
        d = json.load(open(p))
        seeds = d.get("oracle_seeds") or []
        vals = [s["final_de_overlap"]["mean"] for s in seeds
                if s.get("final_de_overlap") is not None]
        if len(vals) < 2:
            continue
        n = len(vals)
        curve = []
        for k in range(1, n + 1):
            curve.append(float(np.mean([max(c) for c in combinations(vals, k)])))
        rows.append({
            "module": d.get("module") or p.stem.replace("step1_", ""),
            "n_seeds": n,
            "seed_values": [round(v, 4) for v in vals],
            "spread": round(max(vals) - min(vals), 4),
            "curve": [round(c, 4) for c in curve],
            # what the last added seed is still buying, in absolute terms and as a
            # fraction of the total climb from one seed to n
            "last_gain": round(curve[-1] - curve[-2], 4),
            "total_climb": round(curve[-1] - curve[0], 4),
        })

    print(f"{'module':<24}{'spread':>8}  " +
          "".join(f"{'k=' + str(k):>8}" for k in range(1, 6)) +
          f"{'last':>9}{'climb':>8}")
    for r in sorted(rows, key=lambda x: -x["last_gain"]):
        c = r["curve"] + [None] * (5 - len(r["curve"]))
        cells = "".join(f"{v:>8.4f}" if v is not None else f"{'-':>8}" for v in c)
        print(f"{r['module']:<24}{r['spread']:>8.4f}  {cells}"
              f"{r['last_gain']:>+9.4f}{r['total_climb']:>8.4f}")

    last = [r["last_gain"] for r in rows]
    climb = [r["total_climb"] for r in rows]
    # the share of the one-to-five climb that the fifth seed alone contributes; a converged
    # search would be spending its last seed on nearly nothing
    share = [r["last_gain"] / r["total_climb"] for r in rows if r["total_climb"] > 0]

    print(f"\nmodules: {len(rows)}")
    print(f"gain from the fifth seed: median {np.median(last):+.4f}, "
          f"max {max(last):+.4f}, still positive on {sum(x > 0 for x in last)} of {len(rows)}")
    print(f"climb from one seed to five: median {np.median(climb):.4f}")
    if share:
        print(f"share of that climb contributed by the fifth seed alone: "
              f"median {np.median(share):.1%}")
    print("\nA converged search spends its last seed on almost nothing. Read the last column "
          "against\nthe climb column: if the fifth seed is still buying a material fraction of "
          "the total,\nthe maximum is not an estimate of the supremum and 'ceiling' is the "
          "wrong word.")

    out = {"n_modules": len(rows), "rows": rows,
           "median_last_gain": float(np.median(last)),
           "median_total_climb": float(np.median(climb)),
           "median_last_share": float(np.median(share)) if share else None,
           "n_still_climbing": int(sum(x > 0 for x in last))}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
