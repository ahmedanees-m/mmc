"""CRISPRi knockdown operator: drive a gene's production to zero, re-simulate."""
from __future__ import annotations
import copy
import numpy as np
from .simulate import steady_state
from ..grammar.model_spec import ModelSpec


def knockdown(spec: ModelSpec, params: dict, gene: str, wt=None) -> np.ndarray:
    """Return Delta = perturbed_ss - wt_ss for all genes."""
    p = copy.deepcopy(params)
    i = spec.genes.index(gene)
    p["basal"] = np.asarray(p["basal"], float).copy()
    p["basal"][i] = 0.0
    p["terms"] = {t: v for t, v in p["terms"].items() if t != gene}  # no production for KD gene
    wt = steady_state(spec, params) if wt is None else wt
    ss = steady_state(spec, p, x0=wt)
    return ss - wt
