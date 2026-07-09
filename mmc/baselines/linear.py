"""Regularized linear baselines (PREREG sections 5 and 7).

Two linear comparators.

persistence is the transfer baseline for Tier A: predict a perturbation's late-state
effect as its measured early-state effect, that is, assume the circuit does not change
across states. A mechanistic model has to beat this to show it has learned the rewiring
rather than copying the training state.

reconstruct is the in-context reconstruction comparator for the Step-5 gate, in the
style of the Zhu low-rank reconstruction (about R 0.39). It represents a perturbation by
its perturbed gene's response signature, how that gene behaves as a target under the
training perturbations, and fits a ridge map from that signature to the full response
vector. A held-out perturbation whose gene was never knocked down is then predicted from
the shared response manifold rather than shrunk to zero.
"""
from __future__ import annotations

import numpy as np


def persistence(train_state_deltas: dict[str, dict[str, float]],
                perts: list[str]) -> dict[str, dict[str, float]]:
    """Predict each perturbation's effect as its measured train-state effect."""
    return {p: dict(train_state_deltas.get(p, {})) for p in perts}


def _signature(deltas: dict[str, dict[str, float]], gene: str,
               sig_perts: list[str]) -> np.ndarray:
    return np.array([deltas.get(q, {}).get(gene, 0.0) for q in sig_perts], dtype=float)


def reconstruct(train_deltas: dict[str, dict[str, float]], genes: list[str],
                train_perts: list[str], held_out_perts: list[str],
                l2: float = 1.0) -> dict[str, dict[str, float]]:
    """Ridge reconstruction of held-out perturbations from response signatures."""
    sig_perts = list(train_perts)
    gi = {g: i for i, g in enumerate(genes)}

    X = np.vstack([_signature(train_deltas, p, sig_perts) for p in train_perts])
    Y = np.zeros((len(train_perts), len(genes)))
    for r, p in enumerate(train_perts):
        for g, v in train_deltas.get(p, {}).items():
            if g in gi:
                Y[r, gi[g]] = v

    A = X.T @ X + l2 * np.eye(X.shape[1])
    W = np.linalg.solve(A, X.T @ Y)                     # (n_signature, n_genes)

    out: dict[str, dict[str, float]] = {}
    for p in held_out_perts:
        pred = _signature(train_deltas, p, sig_perts) @ W
        out[p] = {g: float(pred[gi[g]]) for g in genes}
    return out
