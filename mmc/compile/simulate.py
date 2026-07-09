"""Integrate a compiled model to steady state (WT and perturbed)."""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from .to_ode import build_rhs
from ..grammar.model_spec import ModelSpec


def steady_state(spec: ModelSpec, params: dict, x0=None, T: float = 200.0) -> np.ndarray:
    rhs, _ = build_rhs(spec, params)
    n = len(spec.genes)
    if x0 is None:
        x0 = np.full(n, 0.1)
    sol = solve_ivp(rhs, (0.0, T), x0, method="LSODA", rtol=1e-6, atol=1e-8)
    return sol.y[:, -1]
