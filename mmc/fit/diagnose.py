"""The fit-vs-structure diagnostic gate.

A residual is labelled structural only if no seed recovers its sign or drives it
under tolerance; otherwise it is parametric and withheld, since it is the optimizer's
to resolve, not the structure's. Only structural residuals and the convergence
statistics reach the residual reader, so the reasoning step repairs structure only on
failures a re-fit could not have removed.

Contract: diagnose(spec, observed, n_starts) -> (best_fit,
{(pert, gene): 'structural' | 'parametric'}, convergence_stats).
"""
from __future__ import annotations

import numpy as np

from ..compile.perturb import knockdown
from ..compile.simulate import steady_state
from ..grammar.model_spec import ModelSpec
from .fit_params import _gene_index, multi_fit


def residuals(spec: ModelSpec, params: dict, observed: dict) -> dict:
    """{(perturbation, gene): (predicted_delta, observed_delta)} on the training set."""
    gene_index = _gene_index(spec)
    wt = steady_state(spec, params)
    out: dict[tuple[str, str], tuple[float, float]] = {}
    for pert, deltas in observed.items():
        if pert not in gene_index:
            continue
        d = knockdown(spec, params, pert, wt=wt)
        for gene, obs in deltas.items():
            if gene in gene_index:
                out[(pert, gene)] = (float(d[gene_index[gene]]), float(obs))
    return out


def diagnose(spec: ModelSpec, observed: dict, n_starts: int = 16,
             tol: float = 0.5, **kw) -> tuple[dict, dict, dict]:
    fits = multi_fit(spec, observed, n_starts=n_starts, **kw)
    best = fits[0]
    per_seed = [residuals(spec, f["params"], observed) for f in fits]

    labels: dict[tuple[str, str], str] = {}
    for key in per_seed[0]:
        recoverable = False
        for res in per_seed:
            pred, obs = res[key]
            if (pred > 0) == (obs > 0) or abs(pred - obs) < tol:
                recoverable = True
                break
        labels[key] = "parametric" if recoverable else "structural"

    losses = [f["loss"] for f in fits]
    stats = {
        "n_starts": n_starts,
        "loss_best": float(best["loss"]),
        "loss_median": float(np.median(losses)),
        "loss_spread": float(np.std(losses)),
    }
    return best, labels, stats
