"""Pool the outer folds of the per-fold selection arm and pair at fold level.

The per-module comparison counts how many modules the per-fold arm wins, over a median of
three scoreable outer folds each. That count discards the fold-level pairing the design
provides. This script pools every outer fold across the sound modules and pairs the per-fold
arm against the signal-destroyed arm fold by fold: same module, same fold, same observed
response, two structures. It pairs the same folds against the linear map as well.

A module is sound when its linear map clears the 95th percentile of the random-structure
null, which is the same screen the per-module tables use.

Pooling assumes fold difficulty is comparable across modules. The per-module spread shows it
is not, so the pooled interval is wider than a within-module interval would be.

    python scripts/pool_nested_folds.py --comparator results/comparator/twenty_seeds \
        --offset results/comparator/offset --out results/nested/pooled_folds.json

Reads result JSON only; no model is refit and no data store is needed.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

import numpy as np

N_BOOT = 10_000


def paired_delta(a: np.ndarray, b: np.ndarray, seed: int = 0) -> dict:
    """Mean within-fold difference with a paired bootstrap over folds."""
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[~np.isnan(d)]
    if d.size == 0:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "n": 0, "advantage": False}
    rng = np.random.default_rng(seed)
    boots = d[rng.integers(0, d.size, size=(N_BOOT, d.size))].mean(axis=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"delta": float(d.mean()), "lo": float(lo), "hi": float(hi),
            "n": int(d.size), "advantage": bool(lo > 0)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--comparator", default="results/comparator/twenty_seeds")
    ap.add_argument("--offset", default="results/comparator/offset")
    ap.add_argument("--out", default="results/nested/pooled_folds.json")
    args = ap.parse_args()

    nested_v, permuted_v, linear_v, labels = [], [], [], []
    for p in sorted(glob.glob(f"{args.comparator}/*.json")):
        d = json.load(open(p))
        lin_mean = {r["source"]: r for r in d["table"]}["linear"]["de_overlap_mean"]
        if lin_mean <= d["random_null"]["p95"]:
            continue
        m = d["module"]
        rf = os.path.join(args.offset, f"real_{m}.json")
        pf = os.path.join(args.offset, f"permuted_{m}.json")
        if not (os.path.exists(rf) and os.path.exists(pf)):
            continue

        perts = d["module_summary"]["perts"]
        idx = {g: i for i, g in enumerate(perts)}
        outer = [idx[g] for g in d["oracle_nested"]["outer_perts"] if g in idx]

        real = json.load(open(rf))["per_fold_scores"]
        nested = real["oracle_nested"]
        linear = real["linear"]
        permuted = json.load(open(pf))["per_fold_scores"]["oracle"]

        for i in outer:
            if i >= len(nested) or i >= len(permuted) or i >= len(linear):
                continue
            n, q, lin = nested[i], permuted[i], linear[i]
            if n is None or q is None or n != n or q != q:
                continue
            nested_v.append(n)
            permuted_v.append(q)
            linear_v.append(lin if (lin is not None and lin == lin) else np.nan)
            labels.append(f"{m}:{perts[i]}")

    nested_v = np.array(nested_v, float)
    permuted_v = np.array(permuted_v, float)
    linear_v = np.array(linear_v, float)

    vs_perm = paired_delta(nested_v, permuted_v)
    vs_lin = paired_delta(nested_v, linear_v)

    print(f"pooled outer folds: {len(nested_v)} across "
          f"{len({lab.split(':')[0] for lab in labels})} modules")
    print(f"  per-fold mean {nested_v.mean():.4f}, signal-destroyed mean "
          f"{permuted_v.mean():.4f}, linear mean {np.nanmean(linear_v):.4f}")
    print("\nper-fold minus signal-destroyed, paired at fold level:")
    print(f"  {vs_perm['delta']:+.4f} [{vs_perm['lo']:+.4f}, {vs_perm['hi']:+.4f}]  "
          f"n={vs_perm['n']}  advantage: {vs_perm['advantage']}")
    print("per-fold minus linear, paired at fold level, same folds:")
    print(f"  {vs_lin['delta']:+.4f} [{vs_lin['lo']:+.4f}, {vs_lin['hi']:+.4f}]  "
          f"n={vs_lin['n']}  advantage: {vs_lin['advantage']}")

    out = {"n_folds": int(vs_perm["n"]),
           "n_modules": len({lab.split(":")[0] for lab in labels}),
           "nested_mean": float(nested_v.mean()),
           "permuted_mean": float(permuted_v.mean()),
           "linear_mean": float(np.nanmean(linear_v)),
           "nested_minus_permuted": vs_perm,
           "nested_minus_linear": vs_lin,
           "folds": labels}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
