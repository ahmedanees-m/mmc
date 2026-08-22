"""Amendment A20: the corrected-advantage count re-derived on the seed median.

Section 2.31 deprecated the seed maximum in favour of the median, and recorded that
per-seed per-fold scores had never been stored, which blocked an exact re-derivation. The
comparator now retains them and the 13 modules where the linear arm beats its own random
null have been re-run at five seeds.

The count reported here is the one the argument turns on: on how many of the 13 modules
does the oracle ceiling clear a sound linear baseline, when the ceiling is the seed median
rather than the seed maximum. Both are printed, because the difference between them is the
point.

The oracle is leaky by construction: it selects its structure using the held-out responses
it is then scored on. A module it fails to clear is therefore a strong negative, and a
module it clears is not evidence that any honest procedure could.

    python scripts/step_a20_verdict.py --dir results/a20 --out results/a20_verdict.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import median

import numpy as np

N_BOOT = 10_000


def paired_delta(a, b, seed: int = 0) -> dict:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = ~(np.isnan(a) | np.isnan(b))
    d = a[ok] - b[ok]
    if d.size == 0:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "n": 0, "advantage": False}
    rng = np.random.default_rng(seed)
    boots = d[rng.integers(0, d.size, size=(N_BOOT, d.size))].mean(axis=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"delta": float(d.mean()), "lo": float(lo), "hi": float(hi),
            "n": int(d.size), "advantage": bool(lo > 0)}


def benjamini_hochberg(flags: list[bool]) -> list[bool]:
    """Not applicable to interval decisions; kept out deliberately, see the note below."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/a20")
    ap.add_argument("--out", default="results/a20_verdict.json")
    args = ap.parse_args()

    rows = []
    for p in sorted(Path(args.dir).glob("step1_*.json")):
        d = json.load(open(p))
        seeds = d.get("oracle_seeds") or []
        if not seeds:
            continue
        module = d.get("module") or p.stem.replace("step1_", "")

        per_fold = d.get("per_fold_scores", {})
        linear = np.array([np.nan if v is None else v
                           for v in per_fold.get("linear", [])], float)

        # per-seed per-fold oracle scores, the quantity A20 exists to retain
        mats = []
        for s in seeds:
            v = s.get("per_fold_de_overlap")
            if v is not None:
                mats.append([np.nan if x is None else x for x in v])
        if not mats or linear.size == 0:
            continue
        mat = np.array(mats, float)                     # seeds by folds

        means = [s["final_de_overlap"]["mean"] for s in seeds]
        best_i = int(np.argmax(means))

        # the ceiling as a per-fold vector, under each summary
        oracle_max = mat[best_i]
        oracle_med = np.nanmedian(mat, axis=0)

        null_mean = None
        rn = d.get("random_null") or {}
        for v in (rn.values() if isinstance(rn, dict) else []):
            if isinstance(v, dict) and "mean" in v:
                null_mean = v["mean"]
                break

        rows.append({
            "module": module,
            "n_seeds": len(seeds),
            "n_folds": int(np.sum(~np.isnan(linear))),
            "linear_mean": float(np.nanmean(linear)),
            "oracle_seed_max": float(np.nanmean(oracle_max)),
            "oracle_seed_median": float(np.nanmean(oracle_med)),
            "seed_spread": [round(m, 4) for m in means],
            "vs_linear_max": paired_delta(oracle_max, linear),
            "vs_linear_median": paired_delta(oracle_med, linear),
            "null_mean": null_mean,
        })

    n_max = sum(r["vs_linear_max"]["advantage"] for r in rows)
    n_med = sum(r["vs_linear_median"]["advantage"] for r in rows)

    print(f"{'module':<24}{'folds':>6}{'linear':>8}{'orc.max':>9}{'orc.med':>9}"
          f"{'d(med)':>9}  clears(max/med)")
    for r in sorted(rows, key=lambda x: -x["vs_linear_median"]["delta"]):
        print(f"{r['module']:<24}{r['n_folds']:>6}{r['linear_mean']:>8.4f}"
              f"{r['oracle_seed_max']:>9.4f}{r['oracle_seed_median']:>9.4f}"
              f"{r['vs_linear_median']['delta']:>+9.4f}  "
              f"{str(r['vs_linear_max']['advantage']):<5} / "
              f"{r['vs_linear_median']['advantage']}")

    print(f"\nmodules: {len(rows)} (the sound set, where linear beats its own random null)")
    print(f"oracle clears linear on the seed MAXIMUM: {n_max} of {len(rows)}")
    print(f"oracle clears linear on the seed MEDIAN:  {n_med} of {len(rows)}")
    print("\nNo multiplicity correction is applied to these counts. The decision rule of "
          "section 1.2\nis an interval rule rather than a p-value, and the 13 modules are "
          "the whole sound set\nrather than a screen, so the count is reported as a count.")

    out = {"amendment": "A20", "n_modules": len(rows),
           "n_clears_seed_max": n_max, "n_clears_seed_median": n_med, "rows": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
