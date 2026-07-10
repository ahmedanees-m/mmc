"""Integrate a compiled model to steady state (WT and perturbed)."""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..grammar.model_spec import ModelSpec
from .to_ode import build_rhs


def steady_state(spec: ModelSpec, params: dict, x0=None, T: float = 80.0) -> np.ndarray:
    """Integrate to steady state. Production is bounded (logistic gates) and decay is
    linear, so the system settles on the decay timescale; T well past that with LSODA
    reaches the fixed point. Tolerances are set for a smooth non-stiff right-hand side
    and to keep the fit inner loop affordable."""
    rhs, _ = build_rhs(spec, params)
    n = len(spec.genes)
    if x0 is None:
        x0 = np.full(n, 0.1)
    sol = solve_ivp(rhs, (0.0, T), x0, method="LSODA", rtol=1e-5, atol=1e-7)
    return sol.y[:, -1]
