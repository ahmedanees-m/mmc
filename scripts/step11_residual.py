"""Score every source on the perturbation-specific residual (PREREG_v4 amendment A13).

Section 2.2 established that most of what the linear baseline achieves comes from predicting
the mean training response, and that a structural prediction is the difference between a
clamped and an unclamped fixed point, which cannot express a response shared across
perturbations. The primary metric therefore rewards reproducing a bulk response the model
class forfeits by construction.

This removes that component and rescores. For each fold the shared subspace is the top k
right singular vectors of the training response matrix, built from training perturbations
only, so the subspace never sees the held-out response it decomposes. Both the observation
and every prediction are projected onto the orthogonal complement before scoring, so a
prediction consisting only of the shared program scores at the floor rather than being
compared against a quantity it does not address. The DE set in residual space is the m genes
of largest absolute residual, m being that perturbation's DE count under the existing
threshold, which keeps the metric on its original scale and keeps the tie-break rule of
amendment A2.

Amendment A13 fixed k = 1 as primary with k = 2 and 3 as sensitivity, named the three
possible outcomes in advance, and made this a secondary analysis throughout. The
pre-specified primary metric of section 1.2 is unchanged.

Structures are re-used from a recorded comparator result, so no search is repeated.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mmc.data import module_data, module_extract
from mmc.eval import compare, random_null
from mmc.eval.holdout import topk_overlap
from mmc.grammar.model_spec import ModelSpec

N_NULL = 200


def shared_basis(train_rows: np.ndarray, k: int) -> np.ndarray:
    """Top k right singular vectors of the training responses, uncentred."""
    if train_rows.size == 0:
        return np.zeros((0, train_rows.shape[1]))
    _, _, vt = np.linalg.svd(np.asarray(train_rows, float), full_matrices=False)
    return vt[:max(0, min(k, vt.shape[0]))]


def residualise(v: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Component of v orthogonal to the shared subspace."""
    if basis.shape[0] == 0:
        return np.asarray(v, float)
    v = np.asarray(v, float)
    return v - (v @ basis.T) @ basis


def residual_scores(mod, preds: np.ndarray, k: int, *, seed: int = 0):
    """Per-fold overlap in residual space, plus how much DE mass the projection removed."""
    n = len(mod.perts)
    obs = np.asarray(mod.observed, float)
    out = np.full(n, np.nan)
    removed = np.full(n, np.nan)
    for i in range(n):
        m = int(np.count_nonzero(mod.de_mask[i]))
        if m == 0:
            continue
        train = [j for j in range(n) if j != i]
        basis = shared_basis(obs[train], k)
        r_obs = residualise(obs[i], basis)
        r_pred = residualise(np.asarray(preds[i], float), basis)
        denom = float(np.abs(obs[i]).sum())
        removed[i] = 1.0 - float(np.abs(r_obs).sum()) / denom if denom > 0 else np.nan
        target = set(np.argsort(-np.abs(r_obs))[:m].tolist())
        jac, _ = topk_overlap(r_pred, target, m, seed=seed)
        out[i] = jac
    return out, removed


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", required=True, help="comparator result carrying the specs")
    ap.add_argument("--module-def", default="")
    ap.add_argument("--condition", default="")
    ap.add_argument("--k", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    prior = json.load(open(args.result))
    module = prior["module"]
    condition = args.condition or prior.get("condition", "Stim8hr")

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(module, d["regulators"], d["targets"])
    elif module not in module_extract.MODULES:
        genes = next((s["genes"] for s in (prior.get("specs") or {}).values()
                      if isinstance(s, dict) and s.get("genes")), None)
        if genes is None:
            raise SystemExit(f"{module} is not registered and {args.result} records no "
                             f"gene list to rebuild it from")
        module_extract.register_module(module, genes, genes)

    mod = module_data.build_module_data(module, condition)
    summary = module_data.module_summary(mod)
    print(f"{module} / {condition}: {summary['n_genes']} genes, {summary['n_perts']} "
          f"perturbations, {summary['n_folds_with_de']} scoreable folds")

    preds = {
        "linear": compare.fold_predictions(mod, compare.bind_linear(mod)),
        "mean": compare.fold_predictions(mod, compare.bind_mean(mod)),
        "zero": compare.fold_predictions(mod, compare.bind_zero(mod)),
    }
    for name, spec_json in (prior.get("specs") or {}).items():
        if not isinstance(spec_json, dict) or name.endswith("_nested"):
            continue
        spec = ModelSpec.from_json(json.dumps(spec_json))
        preds[name] = compare.fold_predictions(mod, compare.bind_structural(spec, mod))
    print(f"  sources: {', '.join(sorted(preds))}")

    # the null has to be recomputed in residual space; a metric this different cannot be
    # judged against a null computed on the full response
    rng = np.random.default_rng(0)
    n_edges = len(next(iter((prior.get("specs") or {}).values()), {}).get("edges", [])) or 10
    null_preds = []
    for _ in range(N_NULL):
        spec = random_null.sample_spec(list(mod.genes), list(mod.perts), n_edges, rng)
        null_preds.append(compare.fold_predictions(mod, compare.bind_structural(spec, mod)))
    print(f"  residual null: {N_NULL} random structures at {n_edges} edges")

    out = {"module": module, "condition": condition, "module_summary": summary,
           "n_null": N_NULL, "null_edges": n_edges, "by_k": {}}
    for k in args.k:
        rows = {}
        for name, p in preds.items():
            vals, removed = residual_scores(mod, p, k)
            rows[name] = {"per_fold": vals.tolist(),
                          **compare.bootstrap_mean(vals),
                          "de_mass_removed": float(np.nanmean(removed))}
        null_vals = [float(np.nanmean(residual_scores(mod, p, k)[0])) for p in null_preds]
        ref = np.asarray(rows["linear"]["per_fold"], float)
        for name in rows:
            rows[name]["delta_vs_linear"] = compare.paired_delta(
                np.asarray(rows[name]["per_fold"], float), ref)
        out["by_k"][str(k)] = {
            "sources": rows,
            "null": {"mean": float(np.mean(null_vals)),
                     "p95": float(np.percentile(null_vals, 95)),
                     "max": float(np.max(null_vals))},
        }
        print(f"\nk = {k}, shared component removes "
              f"{rows['linear']['de_mass_removed'] * 100:.1f} percent of DE mass; "
              f"residual null mean {out['by_k'][str(k)]['null']['mean']:.4f}, "
              f"p95 {out['by_k'][str(k)]['null']['p95']:.4f}")
        print(f"  {'source':<18}{'residual':>10}{'vs null':>9}{'delta vs linear':>26}{'adv':>5}")
        for name, r in sorted(rows.items(), key=lambda x: -x[1]["mean"]):
            d = r["delta_vs_linear"]
            print(f"  {name:<18}{r['mean']:>10.4f}"
                  f"{r['mean'] - out['by_k'][str(k)]['null']['mean']:>+9.4f}"
                  f"   {d['delta']:+.4f} [{d['lo']:+.4f}, {d['hi']:+.4f}]"
                  f"{'  YES' if d['advantage'] else '   no':>5}")

    p = Path(args.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(p, "w"), indent=2, default=float)
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
