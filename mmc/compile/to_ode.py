"""Compile a ModelSpec + params into a logistic-additive ODE right-hand side.

dx_i/dt = basal_i + production_i(x) - decay_i * x_i
production_i = sum_k prod_ik * PROD_{j in term_k} logistic(sign_ij * w_ikj * x_j - theta_ik)

Logistic (not Hill) for numerical reliability / to avoid the absorbing off-state
(cf. Belgacem 2026, arXiv:2605.01056).
"""
from __future__ import annotations
import numpy as np
from ..grammar.model_spec import ModelSpec

Params = dict  # {'basal': arr, 'decay': arr, 'terms': {target: [{'prod','w':{reg:val},'theta'}]}}


def logistic(z: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-z))


def build_rhs(spec: ModelSpec, params: Params):
    idx = {g: i for i, g in enumerate(spec.genes)}
    n = len(spec.genes)
    basal = np.asarray(params["basal"], float)
    decay = np.asarray(params["decay"], float)
    term_params = params["terms"]

    def production(x: np.ndarray) -> np.ndarray:
        p = np.zeros(n)
        for tgt, terms in term_params.items():
            i = idx[tgt]
            total = 0.0
            for t in terms:
                gate = 1.0
                for reg, w in t["w"].items():
                    s = spec.edge_sign(reg, tgt) or 1
                    gate *= logistic(s * w * x[idx[reg]] - t["theta"])
                total += t["prod"] * gate
            p[i] = total
        return p

    def rhs(_t, x):
        return basal + production(x) - decay * x

    return rhs, idx
