"""Compile a ModelSpec + params into a logistic-additive ODE right-hand side.

dx_i/dt = basal_i + production_i(x) - decay_i * x_i
production_i = sum_k prod_ik * PROD_{j in term_k} logistic(s_ikj * w_ikj * x_j - theta)

s is the gate sign for regulator j in term k (the term's own sign when set, otherwise
the edge sign). theta is per gate (a per-regulator dict) for full AND, OR, and NOT
logic, or a scalar shared across the term for the monotone additive default. State x
is treated as log-expression, so a predicted knockdown delta compares to an observed
log2 fold change.

Logistic (not Hill) for numerical reliability and to avoid the absorbing off-state
(Belgacem 2026, arXiv:2605.01056).
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
        for tgt, param_terms in term_params.items():
            i = idx[tgt]
            spec_terms = spec.rules[tgt].terms if tgt in spec.rules else []
            total = 0.0
            for st, pt in zip(spec_terms, param_terms):
                gate = 1.0
                theta = pt["theta"]  # scalar (shared) or a per-regulator dict
                for reg, w in pt["w"].items():
                    s = st.signs.get(reg) or spec.edge_sign(reg, tgt) or 1
                    th = theta[reg] if isinstance(theta, dict) else theta
                    gate *= logistic(s * w * x[idx[reg]] - th)
                total += pt["prod"] * gate
            p[i] = total
        return p

    def rhs(_t, x):
        return basal + production(x) - decay * x

    return rhs, idx
