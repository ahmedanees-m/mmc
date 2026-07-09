"""v3 Part-4 precondition: the held-out (leave-one-perturbation-out) prediction gate.

Per Decision_1, the earlier 0.927 on Th2 was in-sample Pearson; the structural backend
fits training data by construction (it recovers known-model data at 1.0), so in-sample
fit is not evidence of prediction, and a discovery harvested from an in-sample fit is
curve-fitting. This runs the user's holdout_eval harness against the structural backend:
leave one perturbation out, refit the structure on the rest, predict the held-out one,
and compare the model with the mean, linear, and zero baselines on ACC_DEG and DE-overlap
with bootstrap CIs. Then edge-ablation-on-held-out (Part-4 necessity) and seed-stability.
With a small perturbation count, DE-overlap and CI separation matter more than a point
estimate on ACC_DEG, which saturates when DE genes move together.

Requires the model client (to obtain the module structure once) and the store.
"""
from __future__ import annotations

import sys

import numpy as np

from mmc.baselines import linear as linear_bl
from mmc.compile import structural
from mmc.data import module_extract
from mmc.eval.holdout import (ModuleData, edge_ablation_holdout, loo_evaluate,
                              print_report, seed_stability)
from mmc.fit import fit_structural
from mmc.loop.run import discover
from mmc.shared import store

FDR = 0.10


def get_spec(module: str, condition: str):
    result = discover(module, condition, backend="structural", max_iters=1, n_starts=6)
    return result.ensemble.best().spec


def build_module(spec, module: str, condition: str) -> ModuleData:
    genes = list(spec.genes)
    gi = {g: i for i, g in enumerate(genes)}
    df = store.module_effects(genes, genes, condition)      # any module gene with KD data
    perts = sorted({p for p in df["perturbation"].tolist() if p in gi})
    pi = {p: i for i, p in enumerate(perts)}
    obs = np.zeros((len(perts), len(genes)))
    fdr = np.ones((len(perts), len(genes)))
    for _, r in df.iterrows():
        p, g = r["perturbation"], r["target_gene"]
        if p in pi and g in gi and p != g:
            obs[pi[p], gi[g]] = float(r["effect_size"])
            f = r["fdr"]
            fdr[pi[p], gi[g]] = float(f) if f == f else 1.0
    return ModuleData(genes, perts, obs, fdr < FDR, spec)


def make_backend():
    def fit_fn(spec, train_perts, train_obs, seed=0):
        observed = {p: {spec.genes[j]: float(train_obs[i, j])
                        for j in range(len(spec.genes)) if spec.genes[j] != p}
                    for i, p in enumerate(train_perts)}
        fits = fit_structural.multi_fit(spec, observed, n_starts=4, max_iter=250)
        return {"spec": spec, "params": fits[0]["params"]}

    def predict_fn(model, pert):
        return np.asarray(structural.knockdown(model["spec"], model["params"], pert))

    return fit_fn, predict_fn


def make_linear(genes):
    def linear_fn(train_perts, train_obs, held):
        td = {p: {genes[j]: float(train_obs[i, j]) for j in range(len(genes))}
              for i, p in enumerate(train_perts)}
        pred = linear_bl.reconstruct(td, genes, list(train_perts), [held])[held]
        return np.array([pred.get(g, 0.0) for g in genes])
    return linear_fn


def main() -> None:
    module = sys.argv[1] if len(sys.argv) > 1 else "Th2_GATA3"
    condition = sys.argv[2] if len(sys.argv) > 2 else "Stim8hr"

    spec = get_spec(module, condition)
    mod = build_module(spec, module, condition)
    print(f"=== held-out prediction gate: {module} in {condition} ===")
    print(f"{len(mod.genes)} genes, {len(mod.perts)} perturbations {mod.perts}, "
          f"{len(spec.edges)} edges, {int(mod.de_mask.sum())} DE entries\n")

    fit_fn, predict_fn = make_backend()
    rep = loo_evaluate(mod, fit_fn, predict_fn, linear_fn=make_linear(mod.genes),
                       n_boot=2000, seed=0)
    print_report(rep)

    print("\nedge ablation on held-out (a required edge lowers held-out ACC_DEG):")
    for r in edge_ablation_holdout(mod, fit_fn, predict_fn, seed=0):
        print(f"  {r['edge']:<18} sign {r['sign']:+d}  full {r['holdout_acc_full']}  "
              f"ablated {r['holdout_acc_ablated']}  drop {r['drop']:+.3f}")

    print("\nseed stability:", seed_stability(mod, fit_fn, predict_fn, seeds=(0, 1, 2)))


if __name__ == "__main__":
    main()
