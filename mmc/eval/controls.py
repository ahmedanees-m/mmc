"""Controls (PREREG section 6).

A negative control that must score near zero and a positive control that must score
high, so a metric that looks good for the wrong reason is caught before it is reported.
"""
from __future__ import annotations

import numpy as np

from .metrics import align, pearson


def shuffled_negative(pred_map: dict[str, dict[str, float]],
                      obs_map: dict[str, dict[str, float]],
                      genes: list[str], seed: int = 0) -> float | None:
    """Score predictions against mismatched perturbations. The pooled correlation should
    sit near zero: a model that scores well here is exploiting a gene-level artefact
    rather than perturbation-specific signal."""
    perts = [p for p in obs_map if p in pred_map]
    if len(perts) < 2:
        return None
    rng = np.random.default_rng(seed)
    shuffled = list(perts)
    while True:
        rng.shuffle(shuffled)
        if all(a != b for a, b in zip(perts, shuffled)):
            break
    P, O = [], []
    for p, q in zip(perts, shuffled):
        pv, _ = align(pred_map[p], {}, genes)
        _, ov = align({}, obs_map[q], genes)
        P.append(pv)
        O.append(ov)
    return pearson(np.concatenate(P), np.concatenate(O))


def wt_positive(obs_map: dict[str, dict[str, float]], genes: list[str]) -> float | None:
    """Score the observed response against itself. Perfect by construction; it confirms
    the metric and the alignment are wired correctly."""
    P, O = [], []
    for p in obs_map:
        v, _ = align(obs_map[p], {}, genes)
        P.append(v)
        O.append(v)
    if not P:
        return None
    return pearson(np.concatenate(P), np.concatenate(O))
