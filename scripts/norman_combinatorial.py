"""Step 4: Norman under the field-standard additive split (PREREG_v4 section 5).

The existing Norman result fit on singles only and predicted doubles. That is the
strictly harder extrapolation, and the mechanistic argument predicts failure there by
construction, so a reviewer who knows the benchmark would note that the hard version
was tested and a null reported. The field's standard split trains on all singles plus
a subset of doubles and predicts the held-out doubles. Both are run; this script is the
standard one.

Two design points follow from the Step 1 result rather than from the original plan.

First, a mean-of-training-doubles baseline is included. On the cytokine module 78
percent of everything the linear map achieved came from predicting the mean training
response, and any protocol that lets a model see some doubles must therefore be
compared against simply averaging them, or a weak model will look strong for capturing
a shared response it did not have to reason about.

Second, the additive baseline is fitted on the training doubles and applied to the
held-out ones, rather than fitted per pair on the pair being predicted. The per-pair
fit used in the compose test sees the answer, which makes it an oracle rather than a
baseline; it is retained under its own name so the two are not confused.

Genetic-interaction subtypes are derived from the fitted coefficients by mechanical
quantile rules, because the published coefficient table is not shipped with the GEO
matrix. The rules follow the definitions recorded in PREREG_norman.md section 2 and are
fixed in PREREG_v4 amendment A7 before this runs.

    python scripts/step4_norman_split.py --folds 4 --out /work/results/step4_norman.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mmc.eval import compare

PB = "/norman/norman_pseudobulk.npz"
DE_T = 0.5
MIN_CELLS = 25
MIN_DE_DOUBLE = 5
TOPK = 50
SEED = 0


def additive_fit(sA, sB, dAB, mask):
    X = np.stack([sA[mask], sB[mask]], axis=1)
    coef, *_ = np.linalg.lstsq(X, dAB[mask], rcond=None)
    pred = coef[0] * sA + coef[1] * sB
    ss_res = float(np.sum((dAB[mask] - pred[mask]) ** 2))
    ss_tot = float(np.sum(dAB[mask] ** 2))
    return pred, (1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0), coef


def de_overlap_at_k(pred, obs, k):
    kk = min(k, pred.size)
    top_pred = set(np.argsort(-np.abs(pred), kind="stable")[:kk].tolist())
    top_obs = set(np.argsort(-np.abs(obs), kind="stable")[:kk].tolist())
    return len(top_pred & top_obs) / kk


def acc_deg(pred, obs, de_idx):
    if len(de_idx) == 0:
        return np.nan
    return float(np.mean(np.sign(pred[de_idx]) == np.sign(obs[de_idx])))


def classify(rows: list[dict]) -> None:
    """Assign a genetic-interaction subtype from the fitted coefficients.

    Norman's definitions, per PREREG_norman.md section 2: synergy is two large
    coefficients, suppression two small, epistasis an asymmetric pair where one single
    accounts for the double and the other contributes almost nothing, and neomorphism a
    large deviation from any additive fit. Thresholds are quantiles of this dataset's
    own coefficient distribution rather than absolute numbers, so they cannot be tuned
    toward an outcome.
    """
    total = np.array([r["c1"] + r["c2"] for r in rows])
    asym = np.array([r["asymmetry"] for r in rows])
    dev = np.array([r["na"] for r in rows])
    lo_total, hi_total = np.quantile(total, [1 / 3, 2 / 3])
    lo_asym = np.quantile(asym, 1 / 3)
    hi_dev = np.quantile(dev, 0.9)

    for r, t, a, d in zip(rows, total, asym, dev):
        if d >= hi_dev:
            r["subtype"] = "neomorphic"
        elif a <= lo_asym and t >= lo_total:
            r["subtype"] = "epistasis"
        elif t <= lo_total:
            r["subtype"] = "suppression"
        elif t >= hi_total:
            r["subtype"] = "synergy"
        else:
            r["subtype"] = "additive"


def load_pairs():
    d = np.load(PB, allow_pickle=True)
    gtype = d["group_type"].astype(str)
    gA, gB = d["gene_A"].astype(str), d["gene_B"].astype(str)
    ncells, logfc = d["n_cells"], d["logfc"].astype(np.float64)

    singles = {gA[i]: logfc[i] for i, t in enumerate(gtype)
               if t == "single" and ncells[i] >= MIN_CELLS}
    rows = []
    for i, t in enumerate(gtype):
        if t != "double" or ncells[i] < MIN_CELLS:
            continue
        if gA[i] not in singles or gB[i] not in singles:
            continue
        sA, sB, dAB = singles[gA[i]], singles[gB[i]], logfc[i]
        de_idx = np.where(np.abs(dAB) >= DE_T)[0]
        if len(de_idx) < MIN_DE_DOUBLE:
            continue
        mask = (np.abs(sA) >= DE_T) | (np.abs(sB) >= DE_T) | (np.abs(dAB) >= DE_T)
        _, r2, coef = additive_fit(sA, sB, dAB, mask)
        c1, c2 = float(coef[0]), float(coef[1])
        hi = max(abs(c1), abs(c2))
        rows.append({
            "A": gA[i], "B": gB[i], "sA": sA, "sB": sB, "dAB": dAB,
            "de_idx": de_idx, "mask": mask,
            "c1": abs(c1), "c2": abs(c2),
            "asymmetry": (min(abs(c1), abs(c2)) / hi) if hi > 0 else 1.0,
            "na": float(np.clip(1.0 - r2, 0.0, 1.0)),
        })
    return singles, rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    singles, rows = load_pairs()
    classify(rows)
    print(f"usable singles {len(singles)}, scorable doubles {len(rows)}")
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["subtype"]] = counts.get(r["subtype"], 0) + 1
    print("subtypes:", counts)

    all_singles = np.stack(list(singles.values()))
    L = 1.2 * float(np.quantile(np.abs(all_singles), 0.999))

    def struct_compose(sA, sB):
        a = np.arctanh(np.clip(sA / L, -0.999, 0.999))
        b = np.arctanh(np.clip(sB / L, -0.999, 0.999))
        return L * np.tanh(a + b)

    rng = np.random.default_rng(SEED)
    order = rng.permutation(len(rows))
    folds = np.array_split(order, args.folds)
    methods = ("structural", "additive_trained", "mean_of_singles", "sum_of_singles",
               "mean_train_doubles", "fitted_additive_oracle", "zero")
    per_pair: dict[str, list] = {m: [None] * len(rows) for m in methods}
    acc_pair: dict[str, list] = {m: [None] * len(rows) for m in methods}

    for f, test_idx in enumerate(folds):
        test = set(int(i) for i in test_idx)
        train = [i for i in range(len(rows)) if i not in test]
        # additive coefficients fitted once on the training doubles, then applied
        X, y = [], []
        for i in train:
            r = rows[i]
            m = r["mask"]
            X.append(np.stack([r["sA"][m], r["sB"][m]], axis=1))
            y.append(r["dAB"][m])
        coef, *_ = np.linalg.lstsq(np.concatenate(X), np.concatenate(y), rcond=None)
        train_mean = np.mean(np.stack([rows[i]["dAB"] for i in train]), axis=0)
        print(f"  fold {f}: train {len(train)}, test {len(test)}, "
              f"coefficients {coef[0]:.3f}, {coef[1]:.3f}")

        for i in test:
            r = rows[i]
            sA, sB, dAB = r["sA"], r["sB"], r["dAB"]
            preds = {
                "structural": struct_compose(sA, sB),
                "additive_trained": coef[0] * sA + coef[1] * sB,
                "mean_of_singles": 0.5 * (sA + sB),
                "sum_of_singles": sA + sB,
                "mean_train_doubles": train_mean,
                "fitted_additive_oracle": additive_fit(sA, sB, dAB, r["mask"])[0],
                "zero": np.zeros_like(dAB),
            }
            for m, p in preds.items():
                per_pair[m][i] = de_overlap_at_k(p, dAB, TOPK)
                acc_pair[m][i] = acc_deg(p, dAB, r["de_idx"])

    out = {
        "protocol": {"folds": args.folds, "topk": TOPK, "de_threshold": DE_T,
                     "saturation_L": L, "n_doubles": len(rows),
                     "n_singles": len(singles)},
        "subtype_counts": counts,
        "subtype_rules": {
            "neomorphic": "additive-fit deviation in the top decile",
            "epistasis": "coefficient asymmetry in the bottom tertile and total "
                         "coefficient at or above the lower tertile",
            "suppression": "total coefficient in the bottom tertile",
            "synergy": "total coefficient in the top tertile",
            "additive": "everything else",
        },
        "overall": {}, "by_subtype": {},
    }

    reference = "additive_trained"
    ref = np.array([v if v is not None else np.nan for v in per_pair[reference]], float)
    for m in methods:
        vals = np.array([v if v is not None else np.nan for v in per_pair[m]], float)
        accs = np.array([v if v is not None else np.nan for v in acc_pair[m]], float)
        out["overall"][m] = {
            "de_overlap": compare.bootstrap_mean(vals),
            "acc_deg": compare.bootstrap_mean(accs),
            "delta_vs_" + reference: compare.paired_delta(vals, ref),
        }

    for subtype in sorted(counts):
        idx = [i for i, r in enumerate(rows) if r["subtype"] == subtype]
        entry = {"n": len(idx)}
        sub_ref = ref[idx]
        for m in methods:
            vals = np.array([per_pair[m][i] for i in idx], float)
            entry[m] = {"de_overlap": compare.bootstrap_mean(vals),
                        "delta_vs_" + reference: compare.paired_delta(vals, sub_ref)}
        out["by_subtype"][subtype] = entry

    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2, default=float)

    print(f"\n{'method':<24}{'DE-overlap@50 [95% CI]':<30}{'delta vs additive_trained'}")
    for m in sorted(methods, key=lambda k: -out["overall"][k]["de_overlap"]["mean"]):
        e = out["overall"][m]
        d = e["delta_vs_" + reference]
        print(f"  {m:<22}{e['de_overlap']['mean']:.4f} "
              f"[{e['de_overlap']['lo']:.4f}, {e['de_overlap']['hi']:.4f}]        "
              f"{d['delta']:+.4f} [{d['lo']:+.4f}, {d['hi']:+.4f}] "
              f"{'YES' if d['advantage'] else 'no'}")

    print("\nprimary hypothesis: structural beats additive_trained on epistasis and "
          "suppression")
    for subtype in ("epistasis", "suppression"):
        e = out["by_subtype"].get(subtype)
        if not e:
            continue
        d = e["structural"]["delta_vs_" + reference]
        print(f"  {subtype:<14} n={e['n']:<4} structural "
              f"{e['structural']['de_overlap']['mean']:.4f}  delta {d['delta']:+.4f} "
              f"[{d['lo']:+.4f}, {d['hi']:+.4f}]  "
              f"{'ADVANTAGE' if d['advantage'] else 'no advantage'}")
    print(f"\nwrote {outpath}")


if __name__ == "__main__":
    main()
