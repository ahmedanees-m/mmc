"""MechPert-style consensus baseline (Tier B only).

MechPert (Novo Nordisk and Oxford) predicts a perturbation's effect from a consensus of
its regulatory neighborhood rather than from the perturbation in isolation. Here, for a
held-out test-state perturbation, predict its effect as a similarity-weighted average of
the visible discovery-subset perturbations' test-state effects. Similarity between two
perturbed genes is the correlation of their response signatures in the training state,
how similarly the two genes behave as targets, so the neighborhood is defined by the
training data and only the discovery effects, which Tier B is allowed to see, enter the
prediction.
"""
from __future__ import annotations

import numpy as np


def _signature(deltas: dict[str, dict[str, float]], gene: str,
               sig_perts: list[str]) -> np.ndarray:
    return np.array([deltas.get(q, {}).get(gene, 0.0) for q in sig_perts], dtype=float)


def _similarity(a: np.ndarray, b: np.ndarray) -> float:
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return 0.0
    return max(float(np.corrcoef(a, b)[0, 1]), 0.0)     # only positive neighbours contribute


def predict(train_deltas: dict[str, dict[str, float]],
            test_visible_deltas: dict[str, dict[str, float]],
            genes: list[str], held_out_perts: list[str], discovery_perts: list[str],
            signature_perts: list[str]) -> dict[str, dict[str, float]]:
    """Consensus prediction for each held-out perturbation from the discovery neighbourhood."""
    out: dict[str, dict[str, float]] = {}
    for h in held_out_perts:
        sh = _signature(train_deltas, h, signature_perts)
        weights, effects = [], []
        for d in discovery_perts:
            weights.append(_similarity(sh, _signature(train_deltas, d, signature_perts)))
            effects.append(test_visible_deltas.get(d, {}))
        total = sum(weights)
        if total > 0:
            out[h] = {g: sum(w * e.get(g, 0.0) for w, e in zip(weights, effects)) / total
                      for g in genes}
        else:
            out[h] = {g: 0.0 for g in genes}
    return out
