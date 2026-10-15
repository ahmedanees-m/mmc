"""Assemble a module's perturbation-response matrix for the evaluation harness.

`module_extract.observed_deltas` returns nested dicts, which is what the fitter
wants. The held-out harness and every Step 1 source want the same content as
aligned arrays with an explicit DE mask, and building that in each script invites
the two paths to drift apart. One builder, used by all of them.

A perturbation is any module gene the atlas knocked down. The on-target self effect
is dropped, since the intervention operator imposes it rather than predicting it.
"""
from __future__ import annotations

import numpy as np

from ..eval.holdout import ModuleData
from ..shared import store
from . import module_extract

DEFAULT_FDR = 0.10


def build_module_data(module: str, condition: str, *, fdr: float = DEFAULT_FDR,
                      spec=None) -> ModuleData:
    """Observed deltas and the DE mask for one module in one condition."""
    genes = module_extract.model_genes(module)
    gi = {g: i for i, g in enumerate(genes)}
    df = store.module_effects(genes, genes, condition)

    perts = sorted({p for p in df["perturbation"].tolist() if p in gi})
    if not perts:
        raise ValueError(f"{module} in {condition} has no perturbed module gene")
    pi = {p: i for i, p in enumerate(perts)}

    obs = np.zeros((len(perts), len(genes)))
    q = np.ones((len(perts), len(genes)))
    for _, row in df.iterrows():
        p, g = row["perturbation"], row["target_gene"]
        if p not in pi or g not in gi or p == g:
            continue
        obs[pi[p], gi[g]] = float(row["effect_size"])
        v = row["fdr"]
        q[pi[p], gi[g]] = float(v) if v == v else 1.0

    return ModuleData(genes, perts, obs, q < fdr, spec)


def module_summary(mod: ModuleData) -> dict:
    """Shape figures reported alongside every result table."""
    per_pert = mod.de_mask.sum(axis=1)
    return {
        "n_genes": len(mod.genes),
        "n_perts": len(mod.perts),
        "n_de_entries": int(mod.de_mask.sum()),
        "de_per_pert_mean": round(float(per_pert.mean()), 3),
        "n_folds_with_de": int((per_pert > 0).sum()),
        "perts": list(mod.perts),
    }
