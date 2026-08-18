"""Step 2's deferred analysis: does the diagnostic predict the ceiling? (PREREG_v4 section 3, section 8)

Section 3 built the identifiability diagnostics and deferred the regression to Step 7,
because four modules cannot support it. This fits the relationship across every module
that has a comparator result, and reports which diagnostic carries the signal.

The quantity being predicted is the ceiling's advantage over the linear baseline, the
paired delta of section 1.2, because that is what a reader wants to know in advance:
given a module, is structural modelling worth attempting on it.

Three things this script refuses to do quietly.

It does not pool the two module sources. Section 8's extraction found that regulon and
co-response modules occupy different parts of the diagnostic range, so source is carried
as a covariate and reported, not averaged away.

It does not report an in-sample fit as if it were predictive. Leave-one-module-out
cross-validation is the headline number; the in-sample R squared is shown beside it so
the gap is visible.

It applies the Benjamini-Hochberg correction of section 1.3 across the module-level
advantage tests, which is the Family A the pre-registration names, and reports how many
survive it rather than how many had a nominal advantage.

    python scripts/step2_regime_fit.py --results 'results/step7/step1_*.json' \
        'results/step1_*_full.json' --out results/step2_regime.json
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

from mmc.eval import compare

DIAGNOSTICS = ("perturbation_specific_ratio", "effective_rank", "leading_pc_fraction",
               "effective_rank_normalised")


def load(paths):
    rows = []
    for pattern in paths:
        for f in sorted(glob.glob(pattern)):
            with open(f) as fh:
                d = json.load(fh)
            if "table" not in d or not d.get("random_null"):
                continue
            t = {r["source"]: r for r in d["table"]}
            if "oracle" not in t or "linear" not in t:
                continue
            dg = d["diagnostics"]
            name = d["module"]
            source = ("regulon" if name.startswith("regulon_")
                      else "coresponse" if name.startswith("coresponse_") else "curated")
            pf = d.get("per_fold_scores", {})
            p = float("nan")
            if "oracle" in pf and "linear" in pf:
                a = np.array([np.nan if v is None else v for v in pf["oracle"]], float)
                b = np.array([np.nan if v is None else v for v in pf["linear"]], float)
                p = compare.permutation_p(a, b)
            rows.append({
                "module": name, "source": source,
                "n_folds": d["module_summary"]["n_folds_with_de"],
                "null": d["random_null"]["mean"],
                "linear": t["linear"]["de_overlap_mean"],
                "oracle": t["oracle"]["de_overlap_mean"],
                "ceiling_advantage": t["oracle"]["delta_vs_linear"],
                "advantage": bool(t["oracle"]["advantage"]),
                "p": p,
                **{k: dg.get(k) for k in DIAGNOSTICS},
            })
    # a module present in both the curated and step7 sets is kept once
    seen = {}
    for r in rows:
        seen[r["module"]] = r
    return list(seen.values())


def fit(x: np.ndarray, y: np.ndarray, g: np.ndarray | None = None) -> dict:
    """Ordinary least squares with an intercept, plus leave-one-out cross-validation.

    Section 8 carries the module source as a covariate rather than pooling it away,
    because regulon and co-response modules sit in different parts of the diagnostic
    range and a pooled slope can be produced by the gap between the two groups alone.
    Passing `g` adds an indicator column for the second level and reports the slope
    holding source fixed.
    """
    ok = ~(np.isnan(x) | np.isnan(y))
    if g is not None:
        g = np.asarray(g)[ok]
    x, y = x[ok], y[ok]
    n = x.size
    if n < 4:
        return {"n": int(n), "r2": float("nan"), "loo_r2": float("nan")}
    def design(xv, gv):
        cols = [np.ones_like(xv), xv]
        if gv is not None and len({*gv.tolist()}) > 1:
            lev = sorted({*gv.tolist()})[1]
            cols.append((gv == lev).astype(float))
        return np.stack(cols, axis=1)

    a = design(x, g)
    coef, *_ = np.linalg.lstsq(a, y, rcond=None)
    pred = a @ coef
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    loo = []
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        ai = design(x[m], None if g is None else g[m])
        ci, *_ = np.linalg.lstsq(ai, y[m], rcond=None)
        xi = design(x[i:i + 1], None if g is None else g[i:i + 1])
        # a fold can drop a source level entirely, leaving the row wider than the fit
        loo.append(float(xi[0, :ci.size] @ ci[:xi.shape[1]]) if xi.shape[1] == ci.size
                   else float(ci[0] + ci[1] * x[i]))
    loo = np.asarray(loo)
    ss_loo = float(((y - loo) ** 2).sum())
    loo_r2 = 1 - ss_loo / ss_tot if ss_tot > 0 else float("nan")

    r = float(np.corrcoef(x, y)[0, 1]) if n > 2 else float("nan")
    return {"n": int(n), "slope": float(coef[1]), "intercept": float(coef[0]),
            "source_term": float(coef[2]) if coef.size > 2 else None,
            "pearson_r": r, "r2": float(r2), "loo_r2": float(loo_r2)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = load(args.results)
    if len(rows) < 4:
        print(f"only {len(rows)} modules with a comparator result; the regression needs "
              f"more before it means anything")
    print(f"{len(rows)} modules")

    # Family A: the module-level test that the ceiling beats linear (section 1.3)
    ps = [r["p"] for r in rows]
    flags = compare.benjamini_hochberg(ps, q=0.05)
    for r, f in zip(rows, flags):
        r["bh_reject"] = bool(f)
    nominal = sum(1 for r in rows if r["advantage"])
    corrected = sum(flags)

    g = np.array([r["source"] for r in rows])
    y = np.array([r["ceiling_advantage"] for r in rows], float)
    fits = {d: fit(np.array([r[d] for r in rows], float), y, g) for d in DIAGNOSTICS}
    # The pre-specified target is the ceiling's margin over linear. On modules where the
    # linear arm sits at or below its own random null that margin carries the baseline's
    # noise rather than a property of the module, so the margin over the null is fitted
    # alongside it. The first remains primary; the second is reported, never substituted.
    y_null = np.array([r["oracle"] - r["null"] for r in rows], float)
    fits_null = {d: fit(np.array([r[d] for r in rows], float), y_null, g)
                 for d in DIAGNOSTICS}
    uninformative = [r["module"] for r in rows if r["linear"] <= r["null"]]

    # Carrying source as a covariate only separates it from the diagnostic if the two
    # sources overlap on that diagnostic. Where they do not, the indicator and the
    # diagnostic are the same column twice and neither coefficient means anything on its
    # own, so the confounding is measured and reported next to every fit.
    collinearity = {}
    if len({*g.tolist()}) > 1:
        lev = sorted({*g.tolist()})[1]
        ind = (g == lev).astype(float)
        for d in DIAGNOSTICS:
            xv = np.array([r[d] for r in rows], float)
            ok = ~np.isnan(xv)
            if ok.sum() < 4 or len({*g[ok].tolist()}) < 2:
                continue
            lo1, hi1 = xv[ok][ind[ok] == 0].min(), xv[ok][ind[ok] == 0].max()
            lo2, hi2 = xv[ok][ind[ok] == 1].min(), xv[ok][ind[ok] == 1].max()
            inter = max(0.0, min(hi1, hi2) - max(lo1, lo2))
            union = max(hi1, hi2) - min(lo1, lo2)
            collinearity[d] = {
                "corr_with_source": float(np.corrcoef(xv[ok], ind[ok])[0, 1]),
                "range_overlap": float(inter / union) if union > 0 else float("nan")}
    best = max((k for k in fits if fits[k].get("loo_r2") == fits[k].get("loo_r2")),
               key=lambda k: fits[k]["loo_r2"], default=None)

    by_source = {}
    for s in sorted({r["source"] for r in rows}):
        sub = [r for r in rows if r["source"] == s]
        by_source[s] = {"n": len(sub),
                        "median_spec": float(np.median([r["perturbation_specific_ratio"]
                                                        for r in sub])),
                        "n_advantage": sum(1 for r in sub if r["advantage"])}

    out = {"n_modules": len(rows), "rows": rows, "fits": fits, "best_diagnostic": best,
           "fits_vs_null": fits_null,
           "linear_uninformative": {"modules": uninformative, "n": len(uninformative)},
           "source_collinearity": collinearity,
           "family_a": {"nominal_advantages": nominal, "bh_rejects": corrected,
                        "q": 0.05, "n_tests": len(rows)},
           "by_source": by_source}
    p = Path(args.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as f:
        json.dump(out, f, indent=2, default=float)

    print(f"\n{'diagnostic':<32}{'n':>4}{'pearson r':>11}{'in-sample R2':>14}{'LOO R2':>9}")
    for d, v in fits.items():
        if v.get("n", 0) < 4:
            print(f"  {d:<30}{v['n']:>4}   too few modules")
            continue
        print(f"  {d:<30}{v['n']:>4}{v['pearson_r']:>11.3f}{v['r2']:>14.3f}{v['loo_r2']:>9.3f}")
    print(f"\nbest by out-of-sample fit: {best}")
    print(f"Family A: {nominal} nominal advantages, {corrected} surviving "
          f"Benjamini-Hochberg at q=0.05 over {len(rows)} tests")
    for s, v in by_source.items():
        print(f"  {s:<12} n={v['n']:<3} median spec/shared {v['median_spec']:.2f}  "
              f"advantages {v['n_advantage']}")
    for d, v in collinearity.items():
        if abs(v["corr_with_source"]) > 0.8 or v["range_overlap"] < 0.2:
            print(f"\n{d}: correlated with source at r={v['corr_with_source']:+.2f} and "
                  f"range overlap {v['range_overlap']:.2f}. The slope and the source term "
                  f"are not separable here and neither should be read alone.")
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
