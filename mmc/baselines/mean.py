"""Mean-of-training-perturbations baseline.

Ahlmann-Eltze, Huber, and Anders (Nature Methods 2025) showed this baseline is hard to
beat: for a held-out perturbation, predict the per-gene mean over the training
perturbations' responses. It carries no perturbation-specific information, so it is the
floor a perturbation-response model has to clear to have earned its complexity.
"""
from __future__ import annotations


def fit(train_deltas: dict[str, dict[str, float]], genes: list[str]) -> dict[str, float]:
    """Per-gene mean response over the training perturbations."""
    profile: dict[str, float] = {}
    for g in genes:
        vals = [d[g] for d in train_deltas.values() if g in d]
        profile[g] = sum(vals) / len(vals) if vals else 0.0
    return profile


def predict(profile: dict[str, float], perts: list[str]) -> dict[str, dict[str, float]]:
    """Predict the same mean profile for every held-out perturbation."""
    return {p: dict(profile) for p in perts}
