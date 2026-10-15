"""Uniform sampling over the bounded grammar, for the Step 1 S3 null (PREREG_v4 section 2.2).

The random ensemble answers a question the other sources cannot: how well does a
structure of this size do on this module by construction, before any proposer has
looked at the data. Every reported source is placed as a percentile of this null,
so a source that scores 0.20 on a module where random structures average 0.19 is
described as what it is.

Sampling is uniform over edge sets subject to the grammar's bounds, with the edge
count matched to whichever source the null is being compared against. Regulators
are drawn from the perturbed genes, because a regulator that is never perturbed
contributes nothing the held-out folds can see.
"""
from __future__ import annotations

import numpy as np

from ..grammar.model_spec import MAX_REGS_PER_TERM, Edge, ModelSpec, Rule, Term


def sample_spec(genes: list[str], perts: list[str], n_edges: int,
                rng: np.random.Generator, *, max_in_degree: int = MAX_REGS_PER_TERM,
                p_repress: float = 0.5) -> ModelSpec:
    """One uniform draw: n_edges signed edges, additive (single-term) rules.

    Single-term rules are the monotone default of the grammar. Logic gates are
    sampled by `sample_spec_gated` when the null needs to cover them; keeping the
    default null additive matches how the proposers actually use the grammar, which
    is overwhelmingly single-term.
    """
    if not perts or not genes:
        raise ValueError("need at least one perturbed gene and one target")
    by_target: dict[str, list[str]] = {}
    edges: list[Edge] = []
    seen: set[tuple[str, str]] = set()
    # Cap the attempt count so a request that the grammar cannot satisfy fails
    # loudly rather than spinning.
    max_attempts = max(1000, n_edges * 200)
    attempts = 0
    while len(edges) < n_edges and attempts < max_attempts:
        attempts += 1
        r = perts[rng.integers(len(perts))]
        t = genes[rng.integers(len(genes))]
        if r == t or (r, t) in seen:
            continue
        if len(by_target.get(t, [])) >= max_in_degree:
            continue
        seen.add((r, t))
        sign = -1 if rng.random() < p_repress else 1
        edges.append(Edge(regulator=r, target=t, sign=sign))
        by_target.setdefault(t, []).append(r)
    if len(edges) < n_edges:
        raise ValueError(
            f"grammar cannot hold {n_edges} edges over {len(genes)} genes at "
            f"in-degree {max_in_degree}; drew {len(edges)}"
        )
    rules = {t: Rule(terms=[Term(regulators=list(rs))]) for t, rs in by_target.items()}
    return ModelSpec(genes=list(genes), edges=edges, rules=rules)


def sample_spec_gated(genes: list[str], perts: list[str], n_edges: int,
                      rng: np.random.Generator, *, p_and: float = 0.3) -> ModelSpec:
    """A draw that also samples the logic, splitting some targets into AND terms.

    Used for the grammar-coverage variant of the null, which checks that the null
    is not artificially low because it never exercises the gates.
    """
    base = sample_spec(genes, perts, n_edges, rng)
    by_target: dict[str, list[Edge]] = {}
    for e in base.edges:
        by_target.setdefault(e.target, []).append(e)
    rules = {}
    for t, es in by_target.items():
        regs = [e.regulator for e in es]
        if len(regs) >= 2 and rng.random() < p_and:
            cut = int(rng.integers(1, len(regs)))
            terms = [Term(regulators=regs[:cut]), Term(regulators=regs[cut:])]
        else:
            terms = [Term(regulators=regs)]
        rules[t] = Rule(terms=terms)
    return ModelSpec(genes=list(base.genes), edges=base.edges, rules=rules)


def _defined(null_values) -> np.ndarray:
    """Drop undefined entries. Accepts None as well as NaN, because the values
    round-trip through JSON where NaN is written as null."""
    return np.asarray([float(x) for x in null_values
                       if x is not None and float(x) == float(x)], float)


def percentile_of(value: float, null_values) -> float:
    """Where a source sits in the null distribution, as a percentile in [0, 100]."""
    v = _defined(null_values)
    if v.size == 0 or value is None or value != value:
        return float("nan")
    return float(100.0 * (v < value).mean())


def summarise_null(null_values) -> dict:
    v = _defined(null_values)
    if v.size == 0:
        return {"n": 0}
    return {
        "n": int(v.size),
        "mean": float(v.mean()),
        "sd": float(v.std(ddof=1)) if v.size > 1 else 0.0,
        "p50": float(np.percentile(v, 50)),
        "p95": float(np.percentile(v, 95)),
        "p99": float(np.percentile(v, 99)),
        "max": float(v.max()),
    }


def _assemble(genes, edges, by_target):
    """Shared tail of every sampler: additive single-term rules over the drawn edges."""
    rules = {t: Rule(terms=[Term(regulators=list(rs))]) for t, rs in by_target.items()}
    return ModelSpec(genes=list(genes), edges=edges, rules=rules)


def _draw(genes, perts, n_edges, rng, pick_regulator, *,
          max_in_degree: int = MAX_REGS_PER_TERM, p_repress: float = 0.5,
          pick_target=None, label: str = "structure"):
    """Draw n_edges under a caller-supplied regulator choice.

    Every topology family in this module differs only in how a regulator is chosen for
    each edge, so the acceptance rules, the in-degree cap, the sign draw and the failure
    behaviour live here rather than being repeated and allowed to drift apart.
    """
    if not perts or not genes:
        raise ValueError("need at least one perturbed gene and one target")
    by_target: dict[str, list[str]] = {}
    edges: list[Edge] = []
    seen: set[tuple[str, str]] = set()
    max_attempts = max(1000, n_edges * 200)
    attempts = 0
    while len(edges) < n_edges and attempts < max_attempts:
        attempts += 1
        r = pick_regulator()
        t = pick_target(r) if pick_target is not None else genes[rng.integers(len(genes))]
        if r == t or (r, t) in seen:
            continue
        if len(by_target.get(t, [])) >= max_in_degree:
            continue
        seen.add((r, t))
        edges.append(Edge(regulator=r, target=t,
                          sign=-1 if rng.random() < p_repress else 1))
        by_target.setdefault(t, []).append(r)
    if len(edges) < n_edges:
        raise ValueError(
            f"grammar cannot hold {n_edges} edges over {len(genes)} genes at "
            f"in-degree {max_in_degree} under {label}; drew {len(edges)}"
        )
    return _assemble(genes, edges, by_target)


def sample_spec_hub(genes: list[str], perts: list[str], n_edges: int,
                    rng: np.random.Generator, *, n_hubs: int = 3,
                    hub_share: float = 0.7, **kw) -> ModelSpec:
    """A few master regulators carry most of the out-edges (PREREG_v4 amendment A14).

    The uniform sampler places regulators evenly, which is not how regulatory networks are
    arranged. This is the direct test of whether a hub-dominated causal generator can
    produce the low-rank signature the uniform ones cannot.
    """
    hubs = [perts[i] for i in rng.permutation(len(perts))[:max(1, min(n_hubs, len(perts)))]]

    def pick():
        if rng.random() < hub_share:
            return hubs[rng.integers(len(hubs))]
        return perts[rng.integers(len(perts))]

    return _draw(genes, perts, n_edges, rng, pick, label="hub", **kw)


def sample_spec_scale_free(genes: list[str], perts: list[str], n_edges: int,
                           rng: np.random.Generator, *, seed_degree: float = 0.2,
                           **kw) -> ModelSpec:
    """Out-degree grown by preferential attachment (PREREG_v4 amendment A14).

    Covers the master-regulator objection without the hub count or their share being
    chosen by hand: a regulator already carrying edges is proportionally more likely to
    take the next one, so the degree distribution is heavy-tailed by construction.
    """
    weight = {p: float(seed_degree) for p in perts}

    def pick():
        tot = sum(weight.values())
        x = rng.random() * tot
        run = 0.0
        for p, w in weight.items():
            run += w
            if x <= run:
                weight[p] += 1.0
                return p
        return perts[-1]

    return _draw(genes, perts, n_edges, rng, pick, label="scale-free", **kw)


def sample_spec_modular(genes: list[str], perts: list[str], n_edges: int,
                        rng: np.random.Generator, *, n_blocks: int = 4,
                        within: float = 0.8, **kw) -> ModelSpec:
    """Edges fall mostly inside gene blocks (PREREG_v4 amendment A14).

    The other common departure from uniformity. Produces block structure rather than a
    dominant direction, so it separates "low rank because one program dominates" from
    "low rank because the module is several weakly coupled programs".
    """
    order = list(rng.permutation(len(genes)))
    block_of = {}
    for rank, gi in enumerate(order):
        block_of[genes[gi]] = rank % max(1, n_blocks)
    by_block: dict[int, list[str]] = {}
    for g, b in block_of.items():
        by_block.setdefault(b, []).append(g)

    def pick():
        return perts[rng.integers(len(perts))]

    def target_for(r):
        b = block_of.get(r)
        if b is not None and rng.random() < within and by_block.get(b):
            pool = by_block[b]
            return pool[rng.integers(len(pool))]
        return genes[rng.integers(len(genes))]

    return _draw(genes, perts, n_edges, rng, pick, pick_target=target_for,
                 label="modular", **kw)


TOPOLOGIES = {
    "uniform": sample_spec,
    "hub": sample_spec_hub,
    "scale_free": sample_spec_scale_free,
    "modular": sample_spec_modular,
}
