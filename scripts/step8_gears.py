"""Step 8: GEARS as a comparator on the combinatorial regime.

The pre-registration scoped this defensively, and the third review made the case for
running it: claim 1 says no structure source beats a linear map, and GEARS is the method a
reader will want that claim tested against, because it injects an external gene-gene graph
that this grammar has no way to represent.

**Comparability, stated before the numbers.** GEARS ships its own processed version of the
Norman atlas, which is not the pseudobulk section 6 was built from. Rather than force one
onto the other, every arm here is computed inside GEARS' own data and split: GEARS itself,
a fitted-additive oracle, and the mean of the two singles. The metric is the project's
DE-overlap, imported from `mmc.eval` rather than reimplemented, so the tie-break rule
cannot drift. The result is internally consistent and is not pooled with section 6's
table, which used a different pseudobulk.

The fitted-additive arm is again the harder bar: it fits its two coefficients against the
observed double it is then scored on, so it sees the answer. GEARS does not.

    python scripts/step8_gears.py --data /data --out /work/results/step8_gears.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mmc.eval.holdout import topk_overlap

TOPK = 50
N_BOOT = 10_000


def de_overlap(pred: np.ndarray, obs: np.ndarray, k: int = TOPK, seed: int = 0) -> float:
    """Precision at k of the predicted top-k against the observed top-k."""
    k = min(k, pred.size)
    if k == 0:
        return float("nan")
    top_o = set(np.argsort(-np.abs(obs), kind="stable")[:k].tolist())
    jac, inter = topk_overlap(pred, top_o, k, seed=seed)
    if jac != jac:
        return float("nan")
    return float(inter / k)


def paired_delta(a: np.ndarray, b: np.ndarray, seed: int = 0) -> dict:
    """PREREG section 1.2: paired bootstrap over folds, advantage only if lo > 0."""
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


def mean_profile(adata, condition: str) -> np.ndarray | None:
    m = adata.obs["condition"].values == condition
    if m.sum() == 0:
        return None
    x = adata.X[m]
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    return np.asarray(x, float).mean(axis=0)


def fitted_additive(sa: np.ndarray, sb: np.ndarray, obs: np.ndarray) -> np.ndarray:
    """Least-squares c1*sa + c2*sb against the observed double. Sees the answer."""
    A = np.vstack([sa, sb]).T
    c, *_ = np.linalg.lstsq(A, obs, rcond=None)
    return A @ c


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/data")
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    from gears import GEARS, PertData

    pert_data = PertData(args.data)
    pert_data.load(data_name="norman")
    pert_data.prepare_split(split="simulation", seed=1)
    pert_data.get_dataloader(batch_size=args.batch_size,
                             test_batch_size=args.batch_size)

    adata = pert_data.adata
    ctrl = mean_profile(adata, "ctrl")
    print(f"cells {adata.shape[0]}, genes {adata.shape[1]}, "
          f"conditions {adata.obs['condition'].nunique()}", flush=True)

    model = GEARS(pert_data, device=args.device)
    model.model_initialize(hidden_size=64)
    print(f"training GEARS for {args.epochs} epochs on {args.device}", flush=True)
    model.train(epochs=args.epochs)

    # test-set doubles whose two singles are both measured, so the additive arms exist
    test = [c for c in pert_data.set2conditions["test"] if "ctrl" not in c.split("+")]
    rows = []
    for cond in test:
        a, b = cond.split("+")
        obs = mean_profile(adata, cond)
        sa = mean_profile(adata, f"{a}+ctrl")
        sb = mean_profile(adata, f"ctrl+{b}")
        if obs is None or sa is None or sb is None:
            continue
        obs_d, sa_d, sb_d = obs - ctrl, sa - ctrl, sb - ctrl

        pred = model.predict([[a, b]])
        key = next(iter(pred))
        gears_d = np.asarray(pred[key], float) - ctrl

        add_d = fitted_additive(sa_d, sb_d, obs_d)
        mean_d = (sa_d + sb_d) / 2.0

        # non-additivity, as in section 6: one minus the additive fit's R squared
        ss_res = float(((obs_d - add_d) ** 2).sum())
        ss_tot = float(((obs_d - obs_d.mean()) ** 2).sum())
        na = 1.0 - (1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0)

        rows.append({
            "condition": cond,
            "non_additivity": na,
            "gears": de_overlap(gears_d, obs_d),
            "fitted_additive": de_overlap(add_d, obs_d),
            "mean_of_singles": de_overlap(mean_d, obs_d),
        })
        print(f"  {cond:<24} na {na:.3f}  gears {rows[-1]['gears']:.3f}  "
              f"fitadd {rows[-1]['fitted_additive']:.3f}  "
              f"mean {rows[-1]['mean_of_singles']:.3f}", flush=True)

    if not rows:
        raise SystemExit("no scoreable test doubles")

    rows.sort(key=lambda r: r["non_additivity"])
    t = max(1, len(rows) // 3)
    sets = {"additive_control": rows[:t], "non_additive": rows[-t:], "all": rows}

    out = {"n_test_doubles": len(rows), "topk": TOPK, "epochs": args.epochs, "sets": {}}
    for name, rs in sets.items():
        g = np.array([r["gears"] for r in rs], float)
        fa = np.array([r["fitted_additive"] for r in rs], float)
        mo = np.array([r["mean_of_singles"] for r in rs], float)
        out["sets"][name] = {
            "n": len(rs),
            "non_additivity_mean": float(np.mean([r["non_additivity"] for r in rs])),
            "gears_mean": float(np.nanmean(g)),
            "fitted_additive_mean": float(np.nanmean(fa)),
            "mean_of_singles_mean": float(np.nanmean(mo)),
            "gears_vs_fitted_additive": paired_delta(g, fa),
            "gears_vs_mean_of_singles": paired_delta(g, mo),
        }
        s = out["sets"][name]
        print(f"\n{name}: n={s['n']}  non-additivity {s['non_additivity_mean']:.3f}")
        print(f"  gears {s['gears_mean']:.4f}  fitted-additive "
              f"{s['fitted_additive_mean']:.4f}  mean-of-singles "
              f"{s['mean_of_singles_mean']:.4f}")
        d = s["gears_vs_fitted_additive"]
        print(f"  gears minus fitted-additive {d['delta']:+.4f} "
              f"[{d['lo']:+.4f}, {d['hi']:+.4f}]  advantage {d['advantage']}")
    out["rows"] = rows

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
