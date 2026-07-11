"""Positive control on the most strongly non-additive Norman pairs (top decile of deviation).

This is the most favorable combinatorial setting for the structural model: the pairs whose
observed double deviates most from an additive model. If the model beats the additive baselines
anywhere held-out, it should be here. Reports the top-decile aggregate (model versus
fitted-additive and mean-of-singles on held-out DE-overlap, with bootstrap CIs) and the count of
individual pairs where the model beats both baselines (a trust-to-correct case).

Deterministic, no model calls. Reuses /norman/norman_pseudobulk.npz. Writes
/norman/norman_positive_control.json.
"""
from __future__ import annotations

import json

import numpy as np

PB = "/norman/norman_pseudobulk.npz"
OUT = "/norman/norman_positive_control.json"
DE_T = 0.5
MIN_CELLS = 25
MIN_DE_DOUBLE = 5
TOPK = 50
N_BOOT = 2000
DECILE = 0.10


def additive_fit(sA, sB, dAB, mask):
    X = np.stack([sA[mask], sB[mask]], axis=1)
    coef, *_ = np.linalg.lstsq(X, dAB[mask], rcond=None)
    pred = coef[0] * sA + coef[1] * sB
    ss_res = float(np.sum((dAB[mask] - pred[mask]) ** 2))
    ss_tot = float(np.sum(dAB[mask] ** 2))
    return pred, (1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0)


def de_overlap(pred, obs, de_idx, k):
    if len(de_idx) == 0:
        return np.nan
    kk = min(k, len(de_idx))
    top_pred = set(np.argsort(-np.abs(pred))[:kk].tolist())
    top_obs = set(np.argsort(-np.abs(obs))[:kk].tolist())
    return len(top_pred & top_obs) / kk


def boot_ci(vals, rng, n_boot=N_BOOT):
    vals = np.asarray([v for v in vals if v == v])
    if len(vals) == 0:
        return (np.nan, np.nan, np.nan)
    means = [np.mean(vals[rng.integers(0, len(vals), len(vals))]) for _ in range(n_boot)]
    return (float(np.mean(vals)), float(np.percentile(means, 2.5)),
            float(np.percentile(means, 97.5)))


def main():
    d = np.load(PB, allow_pickle=True)
    gtype = d["group_type"].astype(str)
    gA = d["gene_A"].astype(str)
    gB = d["gene_B"].astype(str)
    ncells = d["n_cells"]
    logfc = d["logfc"].astype(np.float64)

    singles = {gA[i]: logfc[i] for i, t in enumerate(gtype)
               if t == "single" and ncells[i] >= MIN_CELLS}
    doubles = [(gA[i], gB[i], logfc[i]) for i, t in enumerate(gtype)
               if t == "double" and ncells[i] >= MIN_CELLS
               and gA[i] in singles and gB[i] in singles]

    all_singles = np.stack(list(singles.values()))
    scale = 1.2 * np.quantile(np.abs(all_singles), 0.999)

    def struct_compose(sA, sB):
        a = np.arctanh(np.clip(sA / scale, -0.999, 0.999))
        b = np.arctanh(np.clip(sB / scale, -0.999, 0.999))
        return scale * np.tanh(a + b)

    pairs = []
    for A, B, dAB in doubles:
        sA, sB = singles[A], singles[B]
        de_mask = (np.abs(sA) >= DE_T) | (np.abs(sB) >= DE_T) | (np.abs(dAB) >= DE_T)
        de_double = np.where(np.abs(dAB) >= DE_T)[0]
        if int(de_mask.sum()) < MIN_DE_DOUBLE or len(de_double) < MIN_DE_DOUBLE:
            continue
        p_fit, r2 = additive_fit(sA, sB, dAB, de_mask)
        na = float(np.clip(1.0 - r2, 0.0, 1.0))
        row = {"A": A, "B": B, "na": na,
               "deov_model": de_overlap(struct_compose(sA, sB), dAB, de_double, TOPK),
               "deov_fitadd": de_overlap(p_fit, dAB, de_double, TOPK),
               "deov_mean": de_overlap(0.5 * (sA + sB), dAB, de_double, TOPK)}
        pairs.append(row)

    pairs.sort(key=lambda r: r["na"])
    n = len(pairs)
    k = max(3, round(n * DECILE))
    decile = pairs[-k:]

    rng = np.random.default_rng(0)
    agg = {m: boot_ci([r[f"deov_{m}"] for r in decile], rng)
           for m in ("model", "fitadd", "mean")}
    trust_correct = [f"{r['A']}+{r['B']}" for r in decile
                     if r["deov_model"] > r["deov_fitadd"] and r["deov_model"] > r["deov_mean"]]

    result = {
        "n_pairs_total": n, "decile_n": k,
        "na_range_decile": [round(decile[0]["na"], 3), round(decile[-1]["na"], 3)],
        "deov_model": [round(x, 3) for x in agg["model"]],
        "deov_fitted_additive": [round(x, 3) for x in agg["fitadd"]],
        "deov_mean_of_singles": [round(x, 3) for x in agg["mean"]],
        "trust_correct_pairs": trust_correct,
        "n_trust_correct": len(trust_correct),
        "beats_baselines_aggregate": bool(agg["model"][1] > agg["fitadd"][2]
                                          and agg["model"][1] > agg["mean"][2]),
    }
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print("=== Norman positive control (top decile of non-additivity) ===")
    print(f"decile n={k} of {n}, non-additivity range "
          f"{result['na_range_decile'][0]}-{result['na_range_decile'][1]}")
    print(f"held-out DE-overlap: model {result['deov_model']}  "
          f"fitted-additive {result['deov_fitted_additive']}  "
          f"mean-of-singles {result['deov_mean_of_singles']}")
    print(f"pairs where model beats both baselines (trust-to-correct): "
          f"{result['n_trust_correct']} of {k}  {trust_correct}")
    print(f"model beats both baselines in aggregate (separated CIs): "
          f"{result['beats_baselines_aggregate']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
