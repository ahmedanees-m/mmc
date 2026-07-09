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


# The objective targets what the model is judged on: the differentially expressed genes
# and their direction, not the near-zero bulk that a plain mean-squared error rewards
# (PREREG amendment 2026-07-09). DE gene-perturbation pairs are up-weighted, and a hinge
# penalises predicting the wrong direction on a DE pair.
DE_THRESHOLD = 0.5
DE_WEIGHT = 4.0
SIGN_PENALTY = 0.5


def _loss(spec: ModelSpec, x: np.ndarray, observed: dict, gene_index: dict) -> float:
    params = _unpack(spec, x)
    try:
        wt = steady_state(spec, params)
        if not np.all(np.isfinite(wt)):
            return 1e6
    except Exception:
        return 1e6
    werr, wsum, sign_pen, n_de = 0.0, 0.0, 0.0, 0
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
            if gene not in gene_index:
                continue
            pred = float(d[gene_index[gene]])
            obs = float(obs)
            is_de = abs(obs) >= DE_THRESHOLD
            w = 1.0 + (DE_WEIGHT if is_de else 0.0)
            werr += w * (pred - obs) ** 2
            wsum += w
            if is_de:
                sign_pen += max(0.0, -np.sign(obs) * pred)   # >0 only on a wrong-direction prediction
                n_de += 1
    if wsum == 0:
        return 1e6
    loss = werr / wsum
    if n_de:
        loss += SIGN_PENALTY * (sign_pen / n_de)
    return loss


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


def multi_fit(spec: ModelSpec, observed: dict, n_starts: int = 16,
              n_jobs: int | None = None, **kw) -> list[dict]:
    """Fit from n_starts starts, returned sorted by loss (best first).

    The starts are independent, so they run across processes. n_jobs defaults to one
    per available core, capped at n_starts. Each start is a separate CMA-ES run on a
    small integration, so process-level parallelism turns the multi-start wall time from
    the sum of the starts into roughly one start, which is what makes the full module
    tractable inside the loop.
    """
    import os

    if n_jobs is None:
        n_jobs = min(n_starts, max(1, (os.cpu_count() or 2) - 1))
    if n_jobs <= 1:
        fits = [fit(spec, observed, seed=s, **kw) for s in range(n_starts)]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futures = [ex.submit(fit, spec, observed, seed=s, **kw) for s in range(n_starts)]
            fits = [f.result() for f in futures]
    return sorted(fits, key=lambda f: f["loss"])
