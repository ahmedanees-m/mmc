"""Gradient-based multi-start fitting for the structural steady-state backend (v3 Part 2).

The fit gradient comes from autodiff through the vectorised fixed-point solve (see
compile/structural_jax.py), so L-BFGS-B converges in seconds and the loss is compiled
once and reused across starts. The interface matches fit/diagnose (multi_fit, residuals,
diagnose) so the loop can select this backend without other changes. Multiple starts are
retained for the diagnostic gate, which labels a residual structural only when the best
fit gets a DE gene's direction wrong and most seeds agree.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from ..compile import structural
from ..compile.structural_jax import (bounds, build_tensors, make_loss, param_size,
                                       predict_deltas, to_param_dict)
from ..grammar.model_spec import ModelSpec


def _prep(spec: ModelSpec, observed: dict):
    import jax.numpy as jnp

    t = build_tensors(spec)
    gi = {g: i for i, g in enumerate(spec.genes)}
    perts = [p for p in observed if p in gi]
    obs = np.zeros((len(perts), t["n"]))
    mask = np.zeros((len(perts), t["n"]))
    for pi, p in enumerate(perts):
        for g, v in observed[p].items():
            if g in gi:
                obs[pi, gi[g]] = float(v)
                mask[pi, gi[g]] = 1.0
    pgi = jnp.array([gi[p] for p in perts], dtype=jnp.int32)
    return t, perts, pgi, jnp.asarray(obs), jnp.asarray(mask)


def _fit_prepped(t, pgi, obs, mask, n_starts: int, max_iter: int) -> list[dict]:
    import jax.numpy as jnp

    size = param_size(t)
    if t["T"] == 0 or obs.shape[0] == 0:
        z = np.zeros(size)
        return [{"params": to_param_dict(z, t), "loss": 1e6, "x": z, "seed": 0}]
    vg = make_loss(t, pgi, obs, mask)

    def f(x):
        v, g = vg(jnp.asarray(x))
        return float(v), np.asarray(g, dtype=float)

    lo, hi = bounds(t)
    bnds = list(zip(lo, hi))
    fits = []
    for s in range(n_starts):
        rng = np.random.default_rng(s)
        x0 = lo + rng.random(size) * (hi - lo)
        res = minimize(f, x0, jac=True, method="L-BFGS-B", bounds=bnds,
                       options={"maxiter": max_iter})
        fits.append({"params": to_param_dict(res.x, t), "loss": float(res.fun),
                     "x": np.asarray(res.x), "seed": s})
    return sorted(fits, key=lambda d: d["loss"])


def multi_fit(spec: ModelSpec, observed: dict, n_starts: int = 8,
              max_iter: int = 300, **kw) -> list[dict]:
    t, _perts, pgi, obs, mask = _prep(spec, observed)
    return _fit_prepped(t, pgi, obs, mask, n_starts, max_iter)


def residuals(spec: ModelSpec, params: dict, observed: dict) -> dict:
    """{(perturbation, gene): (predicted, observed)} via the numpy structural model."""
    gi = {g: i for i, g in enumerate(spec.genes)}
    wt = structural.steady_state(spec, params)
    out: dict[tuple[str, str], tuple[float, float]] = {}
    for pert, deltas in observed.items():
        if pert not in gi:
            continue
        d = structural.knockdown(spec, params, pert, wt=wt)
        for gene, obs in deltas.items():
            if gene in gi:
                out[(pert, gene)] = (float(d[gi[gene]]), float(obs))
    return out


def diagnose(spec: ModelSpec, observed: dict, n_starts: int = 8,
             tol: float = 0.5, max_iter: int = 300, **kw) -> tuple[dict, dict, dict]:
    import jax.numpy as jnp

    t, perts, pgi, obs, mask = _prep(spec, observed)
    fits = _fit_prepped(t, pgi, obs, mask, n_starts, max_iter)
    obs_np, mask_np = np.asarray(obs), np.asarray(mask)
    genes = t["genes"]

    per_seed = []
    for fdict in fits:
        deltas = np.asarray(predict_deltas(jnp.asarray(fdict["x"]), t, pgi))
        res = {}
        for pi, p in enumerate(perts):
            for gj, g in enumerate(genes):
                if mask_np[pi, gj] > 0:
                    res[(p, g)] = (float(deltas[pi, gj]), float(obs_np[pi, gj]))
        per_seed.append(res)

    de_tol = 0.5
    n = max(1, len(per_seed))
    labels: dict[tuple[str, str], str] = {}
    for key in (per_seed[0] if per_seed else {}):
        best_pred, o = per_seed[0][key]
        is_de = abs(o) >= de_tol
        best_ok = (best_pred > 0) == (o > 0) or abs(best_pred - o) < tol
        wrong = sum(1 for r in per_seed if (r[key][0] > 0) != (o > 0))
        labels[key] = "structural" if (is_de and not best_ok and wrong >= 0.5 * n) else "parametric"

    best = fits[0]
    losses = [f["loss"] for f in fits]
    stats = {"n_starts": n_starts, "loss_best": float(best["loss"]),
             "loss_median": float(np.median(losses)), "loss_spread": float(np.std(losses))}
    return best, labels, stats
