"""Structural steady-state backend (v3 Part 2).

Model each gene's steady-state level as a bounded structural function of its regulators,
x_i* = g_i(x_pa(i)), with g_i the interpretable logistic sum-of-products already in the
grammar. The steady state is the fixed point x* = F(x*), solved directly by damped
iteration rather than by integrating dynamics, because the data are interventional steady
states and the transients an ODE would add are unconstrained by them. A CRISPRi knockdown
is a do(x_g = 0) intervention: clamp the node to zero and solve for the rest.

This replaces compile/simulate.py + compile/perturb.py for the fit path. The parameters
are basal levels and the per-term production, weight, and threshold; there is no decay or
transient, which makes the fit both faster and better identified.
"""
from __future__ import annotations

import numpy as np

from ..grammar.model_spec import ModelSpec


def _production(spec: ModelSpec, params: dict, x: np.ndarray, idx: dict) -> np.ndarray:
    n = len(spec.genes)
    p = np.zeros(n)
    for tgt, param_terms in params["terms"].items():
        i = idx[tgt]
        spec_terms = spec.rules[tgt].terms if tgt in spec.rules else []
        total = 0.0
        for st, pt in zip(spec_terms, param_terms):
            gate = 1.0
            theta = pt["theta"]
            for reg, w in pt["w"].items():
                s = st.signs.get(reg) or spec.edge_sign(reg, tgt) or 1
                th = theta[reg] if isinstance(theta, dict) else theta
                z = s * w * x[idx[reg]] - th
                gate *= 1.0 / (1.0 + np.exp(-np.clip(z, -40.0, 40.0)))
            total += pt["prod"] * gate
        p[i] = total
    return p


def steady_state(spec: ModelSpec, params: dict, clamp: dict | None = None,
                 iters: int = 400, damp: float = 0.5, tol: float = 1e-7) -> np.ndarray:
    """Solve x* = basal + production(x*) by damped fixed-point iteration. clamp fixes a
    node to a value (a knockdown clamps it to zero)."""
    idx = {g: i for i, g in enumerate(spec.genes)}
    n = len(spec.genes)
    basal = np.asarray(params["basal"], float)
    clamp = clamp or {}
    x = np.full(n, 0.1)
    for g, v in clamp.items():
        x[idx[g]] = v
    for _ in range(iters):
        nxt = basal + _production(spec, params, x, idx)
        for g, v in clamp.items():
            nxt[idx[g]] = v
        nxt = (1.0 - damp) * x + damp * nxt
        if np.max(np.abs(nxt - x)) < tol:
            x = nxt
            break
        x = nxt
    return x


def knockdown(spec: ModelSpec, params: dict, gene: str, wt=None) -> np.ndarray:
    """Delta = perturbed_ss - wt_ss for a do(x_gene = 0) intervention (CRISPRi)."""
    wt = steady_state(spec, params) if wt is None else wt
    ss = steady_state(spec, params, clamp={gene: 0.0})
    return ss - wt


# CRISPRa overexpression clamps the node high rather than to zero. The level is fixed and
# above the bounded WT range (basal is bounded to 3), so it is a genuine gain of function.
ACTIVATION_LEVEL = 4.0


def perturb_set(spec: ModelSpec, params: dict, genes, level: float, wt=None) -> np.ndarray:
    """Delta = perturbed_ss - wt_ss for a do(x_g = level) intervention on a set of genes.
    level 0 is a knockdown; a high level is a CRISPRa overexpression. A set lets a double
    perturbation clamp both genes at once, for the Norman compose test."""
    wt = steady_state(spec, params) if wt is None else wt
    ss = steady_state(spec, params, clamp={g: level for g in genes})
    return ss - wt


def activation(spec: ModelSpec, params: dict, gene: str,
               level: float = ACTIVATION_LEVEL, wt=None) -> np.ndarray:
    """Delta for a single do(x_gene = level) overexpression (CRISPRa)."""
    return perturb_set(spec, params, [gene], level, wt=wt)
