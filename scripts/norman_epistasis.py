"""Work Stream D, steps 2-4: the Norman epistasis compose test (PREREG_norman.md).

Reads the pseudobulk logFC table and runs the pre-registered compose test:

  1. Pair selection (mechanical): for each double (A,B) whose singles are both present,
     score non-additivity as the deviation of the observed double from the best fitted
     additive model c1*sA + c2*sB over the pair's DE genes (1 - R2). Top tertile is the
     non-additive set, bottom tertile the additive control. (PREREG operationalization note:
     the paper's precomputed GI distance-correlation table is not shipped with the GEO
     matrix, so non-additivity is recomputed from the data as the same additive-fit
     deviation. Recorded as a deviation from the letter of PREREG_norman.md section 4.)

  2. Compose test (held-out): a structural logic-gate model (the signed logistic
     sum-of-products kernel, reduced to the bipartite A,B -> readout topology) is informed
     only by the singles and the global response scale, then predicts the double via the
     activation operator (both inputs high), which saturates and so is non-additive. No
     double is seen at fit time. Baselines: fitted-additive (c1*sA + c2*sB fit to the
     observed double, the harder oracle bar) and mean-of-singles.

  3. Gate: DE-overlap (primary) and ACC_DEG on the held-out doubles, bootstrap CIs across
     pairs, per set. The pre-registered win requires the model to beat both additive
     baselines on DE-overlap with separated CIs on the non-additive set and show no such
     advantage on the additive controls.

Writes /norman/norman_result.json and prints a report.
"""
from __future__ import annotations

import json

import numpy as np

PB = "/norman/norman_pseudobulk.npz"
OUT = "/norman/norman_result.json"

DE_T = 0.5          # |logFC| threshold for a differentially expressed gene
MIN_CELLS = 25      # minimum pseudobulk cells for a usable group
MIN_DE_DOUBLE = 5   # a double must move at least this many genes to be scorable
TOPK = 50           # precision-at-k for DE-overlap
N_BOOT = 2000
SEED = 0


def additive_fit(sA, sB, dAB, mask):
    """Least-squares c1*sA + c2*sB ~ dAB over mask; return (pred_full, r2)."""
    X = np.stack([sA[mask], sB[mask]], axis=1)
    y = dAB[mask]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = coef[0] * sA + coef[1] * sB
    ss_res = float(np.sum((y - pred[mask]) ** 2))
    ss_tot = float(np.sum(y ** 2))                     # deltas are relative to control (0)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return pred, r2


def de_overlap(pred, obs, de_idx, k):
    """Precision-at-k: fraction of the top-k predicted DE genes that are truly DE."""
    if len(de_idx) == 0:
        return np.nan
    kk = min(k, len(de_idx))
    top_pred = set(np.argsort(-np.abs(pred))[:kk].tolist())
    top_obs = set(np.argsort(-np.abs(obs))[:kk].tolist())
    return len(top_pred & top_obs) / kk


def acc_deg(pred, obs, de_idx):
    """Sign accuracy over the observed DE genes."""
    if len(de_idx) == 0:
        return np.nan
    return float(np.mean(np.sign(pred[de_idx]) == np.sign(obs[de_idx])))


def boot_ci(vals, rng, n_boot=N_BOOT):
    vals = np.asarray([v for v in vals if v == v])       # drop NaN
    if len(vals) == 0:
        return (np.nan, np.nan, np.nan)
    means = [np.mean(vals[rng.integers(0, len(vals), len(vals))]) for _ in range(n_boot)]
    return (float(np.mean(vals)), float(np.percentile(means, 2.5)),
            float(np.percentile(means, 97.5)))


def main():
    d = np.load(PB, allow_pickle=True)   # self-produced, trusted (string arrays are object dtype)
    genes = d["readout_genes"].astype(str)
    labels = d["group_label"].astype(str)
    gtype = d["group_type"].astype(str)
    gA = d["gene_A"].astype(str)
    gB = d["gene_B"].astype(str)
    ncells = d["n_cells"]
    logfc = d["logfc"].astype(np.float64)                # groups x readout
    n_ro = len(genes)
    print(f"readout genes {n_ro}", flush=True)

    singles, s_ncells = {}, {}
    for i, t in enumerate(gtype):
        if t == "single" and ncells[i] >= MIN_CELLS:
            singles[gA[i]] = logfc[i]
            s_ncells[gA[i]] = int(ncells[i])
    doubles = []
    for i, t in enumerate(gtype):
        if t == "double" and ncells[i] >= MIN_CELLS and gA[i] in singles and gB[i] in singles:
            doubles.append((gA[i], gB[i], logfc[i]))
    print(f"usable singles {len(singles)}  doubles with both singles {len(doubles)}",
          flush=True)

    # global saturation ceiling for the structural (signed-logistic) compose
    all_singles = np.stack(list(singles.values()))
    L = 1.2 * np.quantile(np.abs(all_singles), 0.999)
    print(f"saturation ceiling L = {L:.3f}", flush=True)

    def struct_compose(sA, sB):
        a = np.arctanh(np.clip(sA / L, -0.999, 0.999))
        b = np.arctanh(np.clip(sB / L, -0.999, 0.999))
        return L * np.tanh(a + b)

    pairs = []
    for A, B, dAB in doubles:
        sA, sB = singles[A], singles[B]
        de_mask = (np.abs(sA) >= DE_T) | (np.abs(sB) >= DE_T) | (np.abs(dAB) >= DE_T)
        if int(de_mask.sum()) < MIN_DE_DOUBLE:
            continue
        de_double = np.where(np.abs(dAB) >= DE_T)[0]
        if len(de_double) < MIN_DE_DOUBLE:
            continue
        _, r2 = additive_fit(sA, sB, dAB, de_mask)
        na = float(np.clip(1.0 - r2, 0.0, 1.0))

        p_struct = struct_compose(sA, sB)
        p_mean = 0.5 * (sA + sB)
        p_fit, _ = additive_fit(sA, sB, dAB, de_mask)     # oracle additive (sees the double)

        row = {"A": A, "B": B, "na": na, "n_de": int(len(de_double))}
        for name, pred in (("model", p_struct), ("mean", p_mean), ("fitadd", p_fit)):
            row[f"deov_{name}"] = de_overlap(pred, dAB, de_double, TOPK)
            row[f"acc_{name}"] = acc_deg(pred, dAB, de_double)
        pairs.append(row)

    pairs.sort(key=lambda r: r["na"])
    n = len(pairs)
    t = n // 3
    sets = {"additive_control": pairs[:t], "non_additive": pairs[-t:]}
    print(f"scorable pairs {n}  tertile size {t}", flush=True)

    rng = np.random.default_rng(SEED)
    result = {"n_pairs": n, "tertile_size": t, "L": float(L), "topk": TOPK,
              "de_threshold": DE_T, "sets": {}}
    for sname, rows in sets.items():
        na_vals = [r["na"] for r in rows]
        entry = {"n": len(rows), "na_mean": float(np.mean(na_vals)),
                 "na_range": [float(min(na_vals)), float(max(na_vals))]}
        for metric in ("deov", "acc"):
            entry[metric] = {}
            for method in ("model", "fitadd", "mean"):
                m, lo, hi = boot_ci([r[f"{metric}_{method}"] for r in rows], rng)
                entry[metric][method] = {"mean": m, "ci": [lo, hi]}
        sets[sname] = rows
        result["sets"][sname] = entry

    # pre-registered verdict on DE-overlap (primary)
    def sep(a, b):   # a's CI lower strictly above b's CI upper
        return a["ci"][0] > b["ci"][1]
    na = result["sets"]["non_additive"]["deov"]
    ad = result["sets"]["additive_control"]["deov"]
    win_na = sep(na["model"], na["fitadd"]) and sep(na["model"], na["mean"])
    adv_ad = sep(ad["model"], ad["fitadd"]) and sep(ad["model"], ad["mean"])
    result["verdict"] = {
        "model_beats_both_on_non_additive_separated": bool(win_na),
        "model_advantage_on_additive_controls": bool(adv_ad),
        "prereg_win": bool(win_na and not adv_ad),
    }

    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print("\n=== Norman epistasis compose test ===", flush=True)
    for sname in ("non_additive", "additive_control"):
        e = result["sets"][sname]
        print(f"\n[{sname}] n={e['n']}  non-additivity mean={e['na_mean']:.3f} "
              f"range={e['na_range'][0]:.2f}-{e['na_range'][1]:.2f}")
        for metric, lab in (("deov", "DE-overlap"), ("acc", "ACC_DEG")):
            row = "   ".join(
                f"{m}={e[metric][m]['mean']:.3f}[{e[metric][m]['ci'][0]:.3f},{e[metric][m]['ci'][1]:.3f}]"
                for m in ("model", "fitadd", "mean"))
            print(f"  {lab:11s}: {row}")
    v = result["verdict"]
    print(f"\nPREREG win (model > both additive baselines on non-additive DE-overlap, "
          f"separated CIs, no advantage on additive controls): {v['prereg_win']}")
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
