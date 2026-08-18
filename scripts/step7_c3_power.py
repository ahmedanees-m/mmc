"""Is claim C3 unsupported, or was it untestable with this design?

Section 8's regression returned negative cross-validated R squared on every diagnostic,
which was recorded as no support for C3. That reading is only safe if the design could
have detected a relationship had one been present. Three things argue it might not:

  the target is the ceiling's margin over linear, and linear is at or below its own null
  on 11 of 27 modules, so part of the target is noise;

  the ceiling is a maximum over three search seeds, an upward-biased statistic whose
  variance grows as the seed count falls, and one module spans 0.099 to 0.456 across
  its seeds;

  no power calculation was ever run against the observed scatter at n = 27.

This refits against lower-variance targets and asks what slope the design could have
detected, so the write-up can say which of the two claims it is entitled to make.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np


def cv_r2(x: np.ndarray, y: np.ndarray, g: np.ndarray | None = None) -> float:
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    gg = None if g is None else np.asarray(g)[ok]
    levels = sorted({*gg.tolist()})[1:] if gg is not None else []

    def design(xv, gv):
        cols = [np.ones_like(xv), xv]
        for lev in levels:
            cols.append((np.asarray(gv) == lev).astype(float))
        return np.stack(cols, axis=1)

    n = x.size
    if n < 4:
        return float("nan")
    pred = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        A = design(x[m], None if gg is None else gg[m])
        c, *_ = np.linalg.lstsq(A, y[m], rcond=None)
        x0 = design(x[i:i + 1], None if gg is None else gg[i:i + 1])
        pred[i] = float(x0[0] @ c) if x0.shape[1] == c.size else float(c[0] + c[1] * x[i])
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1 - float(((y - pred) ** 2).sum()) / ss_tot if ss_tot > 0 else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--diagnostic", default="effective_rank_normalised")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-sim", type=int, default=400)
    args = ap.parse_args()

    fit = json.load(open(args.fit))
    rows = {r["module"]: r for r in fit["rows"]}
    seeds, nulls = {}, {}
    for pattern in args.results:
        for f in sorted(glob.glob(pattern)):
            d = json.load(open(f))
            if not (d.get("table") and d.get("random_null")):
                continue
            sp = (d.get("oracle_ceiling_spread") or {}).get("by_seed") or []
            seeds[d["module"]] = [v for v in sp if isinstance(v, (int, float))]
            nulls[d["module"]] = d["random_null"].get("mean")

    names = [m for m in rows if m in seeds and seeds[m]]
    x = np.array([rows[m][args.diagnostic] for m in names], float)
    g = np.array([rows[m]["source"] for m in names])

    targets = {
        "ceiling advantage over linear, max of seeds (as pre-specified)":
            np.array([rows[m]["ceiling_advantage"] for m in names], float),
        "ceiling advantage over linear, mean of seeds":
            np.array([np.mean(seeds[m]) - rows[m]["linear"] for m in names], float),
        "ceiling advantage over linear, median of seeds":
            np.array([np.median(seeds[m]) - rows[m]["linear"] for m in names], float),
        "ceiling over its own null, mean of seeds":
            np.array([np.mean(seeds[m]) - nulls[m] for m in names], float),
        "raw ceiling, mean of seeds":
            np.array([np.mean(seeds[m]) for m in names], float),
    }
    signal = [m for m in names if rows[m]["linear"] > nulls[m]]
    out = {"n": len(names), "diagnostic": args.diagnostic, "targets": {}}
    print(f"{len(names)} modules, diagnostic {args.diagnostic}\n")
    print(f"{'target':<58}{'CV R2':>9}{'resid sd':>10}")
    for label, y in targets.items():
        r2 = cv_r2(x, y, g)
        A = np.stack([np.ones_like(x), x], axis=1)
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        sd = float(np.std(y - A @ c, ddof=2))
        out["targets"][label] = {"cv_r2": r2, "residual_sd": sd, "slope": float(c[1])}
        print(f"  {label:<56}{r2:>9.3f}{sd:>10.4f}")

    # subset where the linear comparator is worth beating at all
    sel = np.array([m in signal for m in names])
    if sel.sum() >= 6:
        y = np.array([np.mean(seeds[m]) - rows[m]["linear"] for m in names], float)
        r2s = cv_r2(x[sel], y[sel], g[sel])
        out["signal_only"] = {"n": int(sel.sum()), "cv_r2": r2s}
        print(f"\nrestricted to the {int(sel.sum())} modules whose linear arm beats its own"
              f" null mean: CV R2 {r2s:.3f}")

    # how large a slope would this design have detected?
    y = targets["ceiling advantage over linear, mean of seeds"]
    A = np.stack([np.ones_like(x), x], axis=1)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    sd = float(np.std(y - A @ c, ddof=2))
    rng = np.random.default_rng(0)
    span = float(x.max() - x.min())
    detect = None
    grid = []
    for mult in np.arange(0.25, 6.01, 0.25):
        slope = mult * sd / span
        hits = 0
        for _ in range(args.n_sim):
            ys = slope * (x - x.mean()) + rng.normal(0, sd, x.size)
            if cv_r2(x, ys, g) > 0:
                hits += 1
        frac = hits / args.n_sim
        grid.append({"slope": float(slope), "rise_over_range": float(slope * span),
                     "power": frac})
        if detect is None and frac >= 0.80:
            detect = grid[-1]
    out["power"] = {"residual_sd": sd, "x_span": span, "grid": grid,
                    "slope_at_80pc_power": detect,
                    "observed_slope": float(c[1]),
                    "observed_rise_over_range": float(c[1] * span)}
    print(f"\npower against the observed scatter, residual sd {sd:.4f}, "
          f"diagnostic range {span:.3f}")
    if detect:
        print(f"  a slope reaches 80 percent power at {detect['slope']:+.4f}, which is a "
              f"rise of {detect['rise_over_range']:.4f} in ceiling advantage")
        print(f"  across the observed range of the diagnostic")
    else:
        print("  no slope in the swept range reaches 80 percent power")
    print(f"  the fitted slope is {c[1]:+.4f}, a rise of {c[1] * span:+.4f} across "
          f"the same range")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=2, default=float)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
