"""Predict the section 8.1 held-out ceilings, before those modules are run.

Section 8.1 requires the five predictions and their intervals to be committed before the
modules are run, so this is deliberately separate from the scoring step and writes a file
that is committed on its own.

Two sets of predictions are produced. The pre-specified one comes from the section 8
regression over all training modules with source carried as a covariate, and supplies the
primary test. The second is fitted on co-response modules alone; three of the five
held-out modules fall outside the co-response training range on two diagnostics, so those
predictions are marked as extrapolated and are reported rather than scored.

The interval is the ordinary least squares prediction interval, which covers the noise in
a new observation as well as the uncertainty in the fitted line, since the criterion asks
where an observed ceiling will fall rather than where its expectation lies.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

DIAGNOSTICS = ("perturbation_specific_ratio", "effective_rank", "leading_pc_fraction",
               "effective_rank_normalised")


def _t_quantile(p: float, df: int) -> float:
    """Student t inverse CDF. SciPy where available, else a standard expansion."""
    try:
        from scipy.stats import t as _t
        return float(_t.ppf(p, df))
    except Exception:
        pass
    # Normal quantile by Acklam's rational approximation, then the Cornish-Fisher
    # correction that carries it to the t distribution. Accurate to about 1e-4 for the
    # degrees of freedom in play here, which is far finer than the interval needs.
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if p < pl:
        q = (-2 * np.log(p)) ** 0.5
        z = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    elif p > ph:
        q = (-2 * np.log(1 - p)) ** 0.5
        z = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    else:
        q = p - 0.5
        r = q * q
        z = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    g1 = (z ** 3 + z) / 4
    g2 = (5 * z ** 5 + 16 * z ** 3 + 3 * z) / 96
    g3 = (3 * z ** 7 + 19 * z ** 5 + 17 * z ** 3 - 15 * z) / 384
    return float(z + g1 / df + g2 / df ** 2 + g3 / df ** 3)


def design(x, g, level):
    cols = [np.ones_like(x), x]
    if g is not None and level is not None:
        cols.append((np.asarray(g) == level).astype(float))
    return np.stack(cols, axis=1)


def fit_and_predict(x, y, g, level, x_new, g_new, conf=0.80):
    A = design(np.asarray(x, float), g, level)
    yv = np.asarray(y, float)
    coef, *_ = np.linalg.lstsq(A, yv, rcond=None)
    resid = yv - A @ coef
    n, p = A.shape
    dof = max(1, n - p)
    s2 = float(resid @ resid) / dof
    XtX_inv = np.linalg.pinv(A.T @ A)
    tq = _t_quantile(0.5 + conf / 2, dof)
    out = []
    for xi, gi in zip(x_new, g_new):
        x0 = design(np.array([float(xi)]), [gi] if g is not None else None, level)[0]
        centre = float(x0 @ coef)
        se = float(np.sqrt(s2 * (1.0 + x0 @ XtX_inv @ x0)))
        out.append({"prediction": centre, "lo": centre - tq * se,
                    "hi": centre + tq * se, "se": se})
    return out, {"n": n, "dof": dof, "coef": [float(c) for c in coef],
                 "residual_sd": float(np.sqrt(s2)), "t_quantile": float(tq)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True, help="step2_regime_fit.py output")
    ap.add_argument("--holdout-diagnostics", required=True)
    ap.add_argument("--target", default="ceiling_advantage",
                    choices=("ceiling_advantage", "ceiling", "ceiling_over_null"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fit = json.load(open(args.fit))
    hold = json.load(open(args.holdout_diagnostics))["modules"]
    rows = fit["rows"]

    # the diagnostic is chosen by cross-validated fit on the training modules only
    chosen = fit.get("best_diagnostic") or DIAGNOSTICS[0]

    def target(r):
        if args.target == "ceiling":
            return r["oracle"]
        if args.target == "ceiling_over_null":
            return r["oracle"] - r["null"]
        return r["ceiling_advantage"]

    x = [r[chosen] for r in rows]
    y = [target(r) for r in rows]
    g = [r["source"] for r in rows]
    level = sorted(set(g))[1] if len(set(g)) > 1 else None

    names = sorted(hold)
    x_new = [hold[m]["diagnostics"][chosen] for m in names]
    g_new = ["coresponse"] * len(names)

    pooled, pooled_meta = fit_and_predict(x, y, g, level, x_new, g_new)

    sub = [(xi, yi) for xi, yi, gi in zip(x, y, g) if gi == "coresponse"]
    within, within_meta = fit_and_predict([a for a, _ in sub], [b for _, b in sub],
                                          None, None, x_new, g_new)
    lo_c = min(a for a, _ in sub); hi_c = max(a for a, _ in sub)

    preds = {}
    for i, m in enumerate(names):
        xi = x_new[i]
        preds[m] = {
            "diagnostic": chosen, "diagnostic_value": xi,
            "pre_specified": pooled[i],
            "within_coresponse": dict(within[i],
                                      extrapolated=bool(xi < lo_c or xi > hi_c)),
        }
    out = {"committed": "before the held-out modules were run",
           "target": args.target, "diagnostic": chosen,
           "training_modules": len(rows),
           "coresponse_training_range": [lo_c, hi_c],
           "pre_specified_fit": pooled_meta, "within_coresponse_fit": within_meta,
           "criterion": "at least 4 of 5 observed values inside the pre-specified "
                        "80 percent interval",
           "predictions": preds}
    p = Path(args.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(p, "w"), indent=2)

    print(f"target {args.target}, diagnostic {chosen}, {len(rows)} training modules")
    print(f"{'module':<22}{'diag':>8}{'predicted':>11}{'80% interval':>24}  within-source")
    for m in names:
        v = preds[m]; ps = v["pre_specified"]; w = v["within_coresponse"]
        print(f"  {m:<20}{v['diagnostic_value']:>8.3f}{ps['prediction']:>11.4f}"
              f"   [{ps['lo']:+.4f}, {ps['hi']:+.4f}]"
              f"   {w['prediction']:+.4f}{'  extrapolated' if w['extrapolated'] else ''}")
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
