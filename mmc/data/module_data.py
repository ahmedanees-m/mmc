"""Assemble a module's perturbation-response matrix for the evaluation harness.

`module_extract.observed_deltas` returns nested dicts, which is what the fitter
wants. The held-out harness and every structure source want the same content as
aligned arrays with an explicit DE mask, and building that in each script invites
the two paths to drift apart. One builder, used by all of them.

A perturbation is any module gene the atlas knocked down. The on-target self effect
is dropped, since the intervention operator imposes it rather than predicting it.

Two input routes. If MMC_MODULE_MATRICES names a directory, the effect-size and FDR
matrices are read from `<module>__<condition>.npz` in it. Otherwise they are queried from
the store named by MMC_ZHU_STORE. Both routes produce the same ModuleData.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from ..eval.holdout import ModuleData
from ..shared import store
from . import module_extract

DEFAULT_FDR = 0.10
MATRICES_ENV = "MMC_MODULE_MATRICES"


def _from_matrices(root: Path, module: str, condition: str, fdr: float,
                   spec) -> ModuleData | None:
    """Read one module's exported matrices, or None when the file is not there."""
    path = root / f"{module}__{condition}.npz"
    if not path.exists():
        return None
    d = np.load(path, allow_pickle=False)
    genes = [str(g) for g in d["genes"]]
    perts = [str(p) for p in d["perts"]]
    return ModuleData(genes, perts, d["effect_size"], d["fdr"] < fdr, spec)


def store_matrices(module: str, condition: str):
    """Query the store for one module: gene labels, perturbation labels, effects, FDRs.

    The FDR matrix is returned unthresholded, so a caller that exports it keeps enough to
    rebuild the DE mask at any threshold.
    """
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

    return genes, perts, obs, q


def build_module_data(module: str, condition: str, *, fdr: float = DEFAULT_FDR,
                      spec=None) -> ModuleData:
    """Observed deltas and the DE mask for one module in one condition."""
    root = os.environ.get(MATRICES_ENV)
    if root:
        from_matrices = _from_matrices(Path(root), module, condition, fdr, spec)
        if from_matrices is not None:
            return from_matrices

    genes, perts, obs, q = store_matrices(module, condition)
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
