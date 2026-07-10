"""JAX structural steady-state engine (v3 Part 2).

The structural model x_i* = g_i(x_pa(i)) is compiled to padded tensors so the fixed point
solves as a vectorised, jit-compiled, damped iteration, and the fit gradient comes from
autodiff through that iteration rather than from finite differences. This is the real
content of the model-class correction: gradient-based fitting on the fixed point, which
is both fast (a fit in seconds) and better identified than integrating ODE dynamics to a
steady state.

A model is (basal, prod, w, theta): a per-gene basal level, and for each product-term a
production scale and a per-regulator-slot weight and threshold. Terms are padded to three
regulator slots and masked. A knockdown is do(x_g = 0): clamp the node and re-solve.
"""
from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np

from ..grammar.model_spec import ModelSpec

jax.config.update("jax_enable_x64", True)
_SLOTS = 3


def build_tensors(spec: ModelSpec) -> dict:
    """Compile a spec to static arrays and the metadata to pack and unpack parameters."""
    gi = {g: i for i, g in enumerate(spec.genes)}
    n = len(spec.genes)
    term_target, term_reg, term_sign, term_mask, term_regnames = [], [], [], [], []
    for tgt, rule in spec.rules.items():
        for term in rule.terms:
            regs = list(term.regulators)[:_SLOTS]
            reg_idx = [gi[r] for r in regs]
            signs = [term.signs.get(r) or spec.edge_sign(r, tgt) or 1 for r in regs]
            pad = _SLOTS - len(regs)
            term_target.append(gi[tgt])
            term_reg.append(reg_idx + [0] * pad)
            term_sign.append(signs + [1] * pad)
            term_mask.append([1.0] * len(regs) + [0.0] * pad)
            term_regnames.append((tgt, regs))
    T = len(term_target)
    return {
        "n": n, "T": T, "genes": list(spec.genes),
        "term_target": jnp.array(term_target, dtype=jnp.int32) if T else jnp.zeros((0,), jnp.int32),
        "term_reg": jnp.array(term_reg, dtype=jnp.int32) if T else jnp.zeros((0, _SLOTS), jnp.int32),
        "term_sign": jnp.array(term_sign, dtype=jnp.float64) if T else jnp.zeros((0, _SLOTS)),
        "term_mask": jnp.array(term_mask, dtype=jnp.float64) if T else jnp.zeros((0, _SLOTS)),
        "regnames": term_regnames,
    }


def param_size(t: dict) -> int:
    return t["n"] + t["T"] + 2 * t["T"] * _SLOTS


def bounds(t: dict) -> tuple[np.ndarray, np.ndarray]:
    n, T = t["n"], t["T"]
    lo = np.concatenate([np.zeros(n), np.zeros(T), np.zeros(T * _SLOTS), np.full(T * _SLOTS, -6.0)])
    hi = np.concatenate([np.full(n, 3.0), np.full(T, 8.0), np.full(T * _SLOTS, 8.0), np.full(T * _SLOTS, 6.0)])
    return lo, hi


def _unpack(vec, t: dict):
    n, T = t["n"], t["T"]
    basal = vec[:n]
    prod = vec[n:n + T]
    w = vec[n + T:n + T + T * _SLOTS].reshape(T, _SLOTS)
    theta = vec[n + T + T * _SLOTS:].reshape(T, _SLOTS)
    return basal, prod, w, theta


def _production(x, prod, w, theta, t):
    xr = x[t["term_reg"]]                              # (T, slots)
    z = t["term_sign"] * w * xr - theta
    gate = jax.nn.sigmoid(z)
    gate = jnp.where(t["term_mask"] > 0, gate, 1.0)
    term_gate = jnp.prod(gate, axis=1)                 # (T,)
    contrib = prod * term_gate
    return jax.ops.segment_sum(contrib, t["term_target"], num_segments=t["n"])


def steady_state(basal, prod, w, theta, t, clamp_mask, clamp_val, iters=100, damp=0.5):
    n = t["n"]
    x0 = jnp.where(clamp_mask > 0, clamp_val, jnp.full(n, 0.1))

    def body(x, _):
        nxt = basal + _production(x, prod, w, theta, t)
        nxt = (1.0 - damp) * x + damp * nxt
        nxt = jnp.where(clamp_mask > 0, clamp_val, nxt)
        return nxt, None

    x, _ = jax.lax.scan(body, x0, None, length=iters)
    return x


def predict_deltas(vec, t, pert_gene_idx, clamp_level: float = 0.0):
    """Predicted deltas (n_perts, n_genes) for do(x_g = clamp_level) at each pert gene.
    clamp_level 0 is a knockdown; a high level is a CRISPRa overexpression (Norman)."""
    basal, prod, w, theta = _unpack(vec, t)
    n = t["n"]
    zeros = jnp.zeros(n)
    wt = steady_state(basal, prod, w, theta, t, zeros, zeros)

    def one(gp):
        cm = jnp.zeros(n).at[gp].set(1.0)
        cv = jnp.zeros(n).at[gp].set(clamp_level)
        ss = steady_state(basal, prod, w, theta, t, cm, cv)
        return ss - wt

    return jax.vmap(one)(pert_gene_idx)


def make_loss(t, pert_gene_idx, obs, obs_mask, de_thr=0.5, de_w=4.0, sign_pen=0.5,
              clamp_level: float = 0.0):
    """Return a jitted value_and_grad of the sign-aware, DE-weighted loss."""
    def loss(vec):
        deltas = predict_deltas(vec, t, pert_gene_idx, clamp_level)
        is_de = (jnp.abs(obs) >= de_thr) & (obs_mask > 0)
        weight = obs_mask * (1.0 + de_w * is_de)
        werr = jnp.sum(weight * (deltas - obs) ** 2) / jnp.maximum(jnp.sum(weight), 1.0)
        wrong = jnp.where(is_de, jnp.maximum(0.0, -jnp.sign(obs) * deltas), 0.0)
        pen = jnp.sum(wrong) / jnp.maximum(jnp.sum(is_de), 1.0)
        return werr + sign_pen * pen

    return jax.jit(jax.value_and_grad(loss))


def to_param_dict(vec, t: dict) -> dict:
    """Convert a fitted flat vector to the {basal, terms} dict the numpy backend uses."""
    vec = np.asarray(vec)
    n, T = t["n"], t["T"]
    basal = np.array(vec[:n], float)
    prod = np.array(vec[n:n + T], float)
    w = np.array(vec[n + T:n + T + T * _SLOTS], float).reshape(T, _SLOTS)
    theta = np.array(vec[n + T + T * _SLOTS:], float).reshape(T, _SLOTS)
    terms: dict[str, list] = {}
    for k, (tgt, regs) in enumerate(t["regnames"]):
        entry = {"prod": float(prod[k]),
                 "w": {r: float(w[k, j]) for j, r in enumerate(regs)},
                 "theta": {r: float(theta[k, j]) for j, r in enumerate(regs)}}
        terms.setdefault(tgt, []).append(entry)
    return {"basal": basal, "terms": terms}
