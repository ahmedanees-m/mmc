"""Summarise the reduced-rank results across modules against two fixed thresholds.

Over the modules where the linear arm beats its own random null:

  primary    the full ridge clears the reduced-rank map on at most four modules
  secondary  the median ratio of reduced-rank to full-rank held-out DE-overlap is at
             least 0.90

Both thresholds are constants below. The script reports whether each holds and does not
adjust them.

    python scripts/summarize_truncation.py --dir results/truncation \
        --out results/truncation/summary.json

Reads the per-module files written by rank_truncation.py; no model is refit and no data
store is needed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import median

CLEARS_ALLOWED = 4          # primary threshold
RATIO_FLOOR = 0.90          # secondary threshold


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/truncation")
    ap.add_argument("--out", default="results/truncation/summary.json")
    args = ap.parse_args()

    rows = []
    for p in sorted(Path(args.dir).glob("*.json")):
        d = json.load(open(p))
        if "tests" not in d:                    # the summary itself, if written in place
            continue
        t = d["tests"]["adaptive"]
        rows.append({
            "module": d["module"],
            "n_perturbations": d["n_perturbations"],
            "effective_rank": d["effective_rank_full"],
            "adaptive_rank": max(set(d["adaptive_ranks"]), key=d["adaptive_ranks"].count),
            "ridge": d["arms"]["full"]["de_overlap"]["mean"],
            "reduced": d["arms"]["adaptive"]["de_overlap"]["mean"],
            "ratio": t["ratio_reduced_over_full"],
            "delta": t["ridge_minus_reduced"]["delta"],
            "lo": t["ridge_minus_reduced"]["lo"],
            "hi": t["ridge_minus_reduced"]["hi"],
            "ridge_clears": t["ridge_clears"],
            "rank_3": d["tests"].get("rank_3", {}).get("ratio_reduced_over_full"),
            "rank_4": d["tests"].get("rank_4", {}).get("ratio_reduced_over_full"),
        })

    n_clears = sum(r["ridge_clears"] for r in rows)
    med = median(r["ratio"] for r in rows)
    primary = n_clears <= CLEARS_ALLOWED
    secondary = med >= RATIO_FLOOR

    print(f"{'module':<24}{'n':>4}{'effR':>7}{'k':>4}{'ridge':>8}{'reduced':>9}"
          f"{'ratio':>7}{'delta':>9}  clears")
    for r in sorted(rows, key=lambda x: x["ratio"]):
        print(f"{r['module']:<24}{r['n_perturbations']:>4}{r['effective_rank']:>7.2f}"
              f"{r['adaptive_rank']:>4}{r['ridge']:>8.4f}{r['reduced']:>9.4f}"
              f"{r['ratio']:>7.3f}{r['delta']:>+9.4f}  {r['ridge_clears']}")

    print(f"\nmodules: {len(rows)}")
    print(f"ridge clears the reduced-rank map on {n_clears}, threshold at most "
          f"{CLEARS_ALLOWED}: primary {'holds' if primary else 'fails'}")
    print(f"median ratio {med:.3f}, threshold at least {RATIO_FLOOR}: "
          f"secondary {'holds' if secondary else 'fails'}")

    out = {"n_modules": len(rows), "rows": rows,
           "rule": {"clears_allowed": CLEARS_ALLOWED, "ratio_floor": RATIO_FLOOR},
           "n_ridge_clears": n_clears, "median_ratio": med,
           "primary_holds": primary, "secondary_holds": secondary}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
