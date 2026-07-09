"""CP-A2 + CP-A3: run the discovery loop on the powered cytokine module and read the
held-out prediction gate.

Run after scripts/cytokine_module.py has written cytokine_module_def.json. The loop
proposes regulator-to-cytokine structure, fits it with the structural backend, reads the
residuals, and repairs; the frozen structure is then scored leave-one-perturbation-out
against the mean, linear, and zero baselines on ACC_DEG and DE-overlap. The success rule
(PREREG_discovery.md) is model greater than linear on held-out DE-overlap with separated
CIs. Dark-candidate edges in the frozen structure are the discovery candidates for CP-A4.
Requires the model client and the store.
"""
from __future__ import annotations

import json

import numpy as np

from mmc.baselines import linear as linear_bl
from mmc.compile import structural
from mmc.data import module_extract
from mmc.eval.holdout import ModuleData, loo_evaluate, print_report, seed_stability
from mmc.fit import fit_structural
from mmc.loop import run as runmod
from mmc.loop.propose import model_hash
from mmc.loop.run import discover
from mmc.shared import store

MODULE, COND, FDR = "Cytokine_production", "Stim8hr", 0.10
CONTEXT = (
    "These genes are candidate regulators of cytokine production in stimulated primary "
    "human CD4+ T cells together with the cytokine and chemokine outputs (IL2, IL4, IL5, "
    "IL13, IFNG, IL17A, IL17F, IL21, IL10, TNF, CSF2 and related). Propose signed "
    "regulator-to-cytokine edges, and a small number of regulator-to-regulator scaffold "
    "edges where one regulator controls another. Knockdown effects are log2 fold changes "
    "at 8 hours of stimulation."
)


def build_module_data(spec):
    genes = list(spec.genes)
    gi = {g: i for i, g in enumerate(genes)}
    df = store.module_effects(genes, genes, COND)
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


def main():
    with open("/app/cytokine_module_def.json") as f:
        d = json.load(f)
    module_extract.register_module(MODULE, d["regulators"], d["targets"])
    runmod.DEFAULT_CONTEXT[MODULE] = CONTEXT
    dark = set(d["dark"])

    print("=== CP-A2: discovery loop on the cytokine module ===")
    result = discover(MODULE, COND, backend="structural", max_iters=4, n_starts=6, max_iter=250)
    frozen = result.ensemble.best()
    for h in result.history:
        print(f"  fit {h['hash']} loss {h['loss']:.4f} structural residuals {h['n_structural']}")
    print("engineer trace:")
    for s in result.trace.steps:
        print(f"  [{s.kind}] {s.model_hash} {s.rationale}")
    dark_edges = [(e.regulator, e.target, e.sign) for e in frozen.spec.edges if e.regulator in dark]
    print(f"frozen {model_hash(frozen.spec)}: {len(frozen.spec.edges)} edges, "
          f"{len(dark_edges)} from dark candidates:")
    for r, t, s in dark_edges:
        print(f"    {r} -> {t} ({s:+d})")

    print("\n=== CP-A3: held-out prediction gate (the hinge) ===")
    mod = build_module_data(frozen.spec)
    print(f"{len(mod.genes)} genes, {len(mod.perts)} perturbations, {int(mod.de_mask.sum())} DE entries")
    fit_fn, predict_fn = make_backend()
    rep = loo_evaluate(mod, fit_fn, predict_fn, linear_fn=make_linear(mod.genes), n_boot=2000, seed=0)
    print_report(rep)
    print("seed stability:", seed_stability(mod, fit_fn, predict_fn, seeds=(0, 1, 2)))


if __name__ == "__main__":
    main()
