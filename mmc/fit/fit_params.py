"""Multi-start parameter fitting (CMA-ES).

Given the structure (a ModelSpec), fit the continuous parameters to the training
knockdown responses only, from several starts, with bounded parameters for stable
integration. The multiple starts let the diagnostic gate (diagnose.py) separate a fit
failure from a structure failure. The optimizer searches over per-gene basal and
decay and, for each rule term, its production, per-regulator weight, and per-regulator
threshold, minimizing the mean squared error between the predicted and observed
knockdown delta on the training set. It returns every seed for diagnose().

The observed deltas are log2 fold changes and the model state is log-expression, so
the log-base constant is absorbed into the fitted parameters.
"""
from __future__ import annotations

import numpy as np

from ..compile.perturb import knockdown
from ..compile.simulate import steady_state
from ..grammar.model_spec import ModelSpec


def _bounds(spec: ModelSpec) -> tuple[np.ndarray, np.ndarray]:
    n = len(spec.genes)
    lo: list[float] = [0.0] * n + [0.2] * n          # basal, decay
    hi: list[float] = [3.0] * n + [3.0] * n
    for _tgt, rule in spec.rules.items():
        for term in rule.terms:
            lo.append(0.0); hi.append(8.0)           # prod
            for _reg in term.regulators:
                lo += [0.0, -6.0]; hi += [8.0, 6.0]  # weight, threshold
    return np.array(lo), np.array(hi)


def _unpack(spec: ModelSpec, x: np.ndarray) -> dict:
    i, n = 0, len(spec.genes)
    basal = np.array(x[i:i + n], float); i += n
    decay = np.array(x[i:i + n], float); i += n
    terms: dict[str, list] = {}
    for tgt, rule in spec.rules.items():
        tlist = []
        for term in rule.terms:
            prod = float(x[i]); i += 1
            w: dict[str, float] = {}
            theta: dict[str, float] = {}
            for reg in term.regulators:
                w[reg] = float(x[i]); i += 1
                theta[reg] = float(x[i]); i += 1
            tlist.append({"prod": prod, "w": w, "theta": theta})
        terms[tgt] = tlist
    return {"basal": basal, "decay": decay, "terms": terms}


def _gene_index(spec: ModelSpec) -> dict[str, int]:
    return {g: i for i, g in enumerate(spec.genes)}


def _loss(spec: ModelSpec, x: np.ndarray, observed: dict, gene_index: dict) -> float:
    params = _unpack(spec, x)
    try:
        wt = steady_state(spec, params)
        if not np.all(np.isfinite(wt)):
            return 1e6
    except Exception:
        return 1e6
    total, count = 0.0, 0
    for pert, deltas in observed.items():
        if pert not in gene_index:
            continue
        try:
            d = knockdown(spec, params, pert, wt=wt)
        except Exception:
            return 1e6
        if not np.all(np.isfinite(d)):
            return 1e6
        for gene, obs in deltas.items():
            if gene in gene_index:
                total += (float(d[gene_index[gene]]) - float(obs)) ** 2
                count += 1
    return total / count if count else 1e6


def fit(spec: ModelSpec, observed: dict, seed: int = 0,
        max_iter: int = 80, sigma0: float = 0.3) -> dict:
    """Fit parameters from one start. Returns {params, loss, x, seed}."""
    import cma

    lo, hi = _bounds(spec)
    gene_index = _gene_index(spec)
    rng = np.random.default_rng(seed)
    x0 = lo + rng.random(len(lo)) * (hi - lo)
    es = cma.CMAEvolutionStrategy(list(x0), sigma0, {
        "bounds": [list(lo), list(hi)], "maxiter": max_iter,
        "verbose": -9, "seed": seed + 1,
    })
    while not es.stop():
        sols = es.ask()
        es.tell(sols, [_loss(spec, np.asarray(s), observed, gene_index) for s in sols])
    xbest = np.asarray(es.result.xbest)
    return {"params": _unpack(spec, xbest), "loss": float(es.result.fbest),
            "x": xbest, "seed": seed}


def multi_fit(spec: ModelSpec, observed: dict, n_starts: int = 16, **kw) -> list[dict]:
    """Fit from n_starts starts, returned sorted by loss (best first)."""
    fits = [fit(spec, observed, seed=s, **kw) for s in range(n_starts)]
    return sorted(fits, key=lambda f: f["loss"])
