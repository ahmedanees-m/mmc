"""The Step 1 S5 ceiling estimator (PREREG_v4 section 2.2).

This searches the grammar for the structure with the best held-out score, using the
held-out score itself as the search objective. That is leaky on purpose. The
question S5 answers is not "what would an honest procedure achieve" but "what is
the best any structure in this grammar could possibly do on this data, even with
access to the answer". If that upper bound does not clear the linear baseline, no
proposer can, and the ceiling is a property of the data.

Because it is leaky, the number it produces is labelled a ceiling estimator
everywhere it appears, and `nested_search` is run alongside it: structure selected
on an inner split, scored on an outer split the selection never touched. The pair
brackets the truth, the leaky number from above and the nested number from below.

The search is greedy forward-backward over a pool of candidate edges ranked by
marginal association, followed by simulated annealing whose moves are drawn from
the full edge space rather than the pool. The annealing phase is what keeps the
pool from silently capping the ceiling: if the pool were the binding constraint,
annealing would find edges outside it and the reported improvement would show that.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..grammar.model_spec import MAX_REGS_PER_TERM, Edge, ModelSpec, Rule, Term

EdgeKey = tuple[str, str, int]     # (regulator, target, sign)


@dataclass
class SearchTrace:
    """Every structure the search scored, for the Step 2 equivalence-width diagnostic."""

    edge_sets: list[tuple[EdgeKey, ...]] = field(default_factory=list)
    holdout: list[float] = field(default_factory=list)
    train_loss: list[float] = field(default_factory=list)

    def add(self, edges, holdout: float, train_loss: float) -> None:
        self.edge_sets.append(tuple(sorted(edges)))
        self.holdout.append(float(holdout))
        self.train_loss.append(float(train_loss))

    def best_index(self) -> int:
        return int(np.nanargmax(self.holdout))

    def to_dict(self) -> dict:
        return {"n_evaluated": len(self.holdout),
                "holdout": [round(v, 4) for v in self.holdout],
                "train_loss": [round(v, 5) for v in self.train_loss],
                "edge_sets": [[f"{r}->{t}:{s:+d}" for r, t, s in es]
                              for es in self.edge_sets]}


def spec_from_edges(genes: list[str], edges) -> ModelSpec:
    """Build an additive (single-term) spec from an edge set."""
    by_target: dict[str, list[str]] = {}
    kept: list[Edge] = []
    for reg, tgt, sign in sorted(edges):
        if len(by_target.get(tgt, [])) >= MAX_REGS_PER_TERM:
            continue
        kept.append(Edge(regulator=reg, target=tgt, sign=int(sign)))
        by_target.setdefault(tgt, []).append(reg)
    rules = {t: Rule(terms=[Term(regulators=rs)]) for t, rs in by_target.items()}
    return ModelSpec(genes=list(genes), edges=kept, rules=rules)


def candidate_pool(mod, size: int = 150) -> list[tuple[str, str]]:
    """Rank possible (regulator, target) pairs by marginal association.

    A pair whose knockdown never moves the target cannot help a held-out
    prediction, so ranking by the observed absolute effect concentrates the greedy
    phase where the data has something to say. The annealing phase is not
    restricted to this pool.
    """
    genes, perts = list(mod.genes), list(mod.perts)
    gi = {g: i for i, g in enumerate(genes)}
    scored = []
    for pi, p in enumerate(perts):
        for g in genes:
            if g == p:
                continue
            eff = abs(float(mod.observed[pi, gi[g]]))
            de = bool(mod.de_mask[pi, gi[g]])
            scored.append(((p, g), eff + (1.0 if de else 0.0)))
    scored.sort(key=lambda kv: -kv[1])
    return [pair for pair, _ in scored[:size]]


def all_pairs(mod) -> list[tuple[str, str]]:
    return [(p, g) for p in mod.perts for g in mod.genes if p != g]


def _in_degree(edges) -> dict[str, int]:
    d: dict[str, int] = {}
    for _, t, _s in edges:
        d[t] = d.get(t, 0) + 1
    return d


def _legal_add(edges, pair, sign) -> bool:
    reg, tgt = pair
    if any(e[0] == reg and e[1] == tgt for e in edges):
        return False
    return _in_degree(edges).get(tgt, 0) < MAX_REGS_PER_TERM


def _rank_and_score(sets, score_many, screen_many, screen_keep, trace):
    """Score a batch of candidates, optionally screening cheaply first.

    Scoring every candidate at full fold count is what makes the search expensive:
    a greedy step over a 120-pair pool is 240 candidates, and at the module sizes
    here that is minutes per step. Screening ranks the batch with a cheaper fold
    count and sends only the top `screen_keep` through the real objective. The
    accepted score is always the real objective, so screening changes which
    candidates get considered, never what a candidate is worth.
    """
    if screen_many is None or len(sets) <= screen_keep:
        scored = score_many(sets)
        if trace is not None:
            for es, (h, ln) in zip(sets, scored):
                trace.add(es, h, ln)
        return sets, scored
    rough = screen_many(sets)
    order = np.argsort([-(h if not np.isnan(h) else -np.inf) for h, _ in rough])
    short = [sets[i] for i in order[:screen_keep]]
    scored = score_many(short)
    if trace is not None:
        for es, (h, ln) in zip(short, scored):
            trace.add(es, h, ln)
    return short, scored


def greedy_forward_backward(score_many, pool_pairs, *, max_edges: int = 30,
                            min_gain: float = 1e-4, trace: SearchTrace | None = None,
                            init: set[EdgeKey] | None = None,
                            screen_many=None, screen_keep: int = 20) -> set[EdgeKey]:
    """Add the best edge until nothing helps, then drop any edge whose removal helps.

    score_many(list_of_edge_sets) returns a list of (holdout, train_loss) pairs, so
    the caller decides whether the candidates at each step are scored serially or
    across a process pool. screen_many has the same shape and is used only to
    shortlist; see `_rank_and_score`.
    """
    current: set[EdgeKey] = set(init or set())
    cur_score, cur_loss = score_many([frozenset(current)])[0]
    if trace is not None:
        trace.add(current, cur_score, cur_loss)

    improved = True
    while improved:
        improved = False

        # forward
        while len(current) < max_edges:
            cands = [(pair, s) for pair in pool_pairs for s in (1, -1)
                     if _legal_add(current, pair, s)]
            if not cands:
                break
            sets = [frozenset(current | {(p[0], p[1], s)}) for p, s in cands]
            sets, scored = _rank_and_score(sets, score_many, screen_many, screen_keep, trace)
            k = int(np.nanargmax([h for h, _ in scored]))
            best_h, best_l = scored[k]
            if not (best_h > cur_score + min_gain):
                break
            current = set(sets[k])
            cur_score, cur_loss = best_h, best_l
            improved = True

        # backward
        if current:
            sets = [frozenset(current - {e}) for e in list(current)]
            sets, scored = _rank_and_score(sets, score_many, screen_many, screen_keep, trace)
            k = int(np.nanargmax([h for h, _ in scored]))
            if scored[k][0] > cur_score + min_gain:
                current = set(sets[k])
                cur_score, cur_loss = scored[k]
                improved = True

    return current


def simulated_annealing(score_many, start: set[EdgeKey], full_pairs, *,
                        n_steps: int = 1500, t0: float = 0.02, t1: float = 0.001,
                        batch: int = 32, max_edges: int = 30, seed: int = 0,
                        trace: SearchTrace | None = None) -> set[EdgeKey]:
    """Metropolis over add / drop / flip-sign moves drawn from the full edge space.

    Proposals are evaluated in batches so the pool stays busy; within a batch the
    best accepted proposal is taken, which is a greedy-within-batch variant of
    Metropolis and behaves the same at the temperatures used here.
    """
    rng = np.random.default_rng(seed)
    current = set(start)
    cur_score = score_many([frozenset(current)])[0][0]
    best, best_score = set(current), cur_score

    n_batches = max(1, n_steps // batch)
    for b in range(n_batches):
        temp = t0 * (t1 / t0) ** (b / max(1, n_batches - 1))
        proposals = []
        for _ in range(batch):
            move = rng.random()
            cand = set(current)
            if cand and move < 0.35:
                cand.discard(list(cand)[rng.integers(len(cand))])
            elif cand and move < 0.60:
                e = list(cand)[rng.integers(len(cand))]
                cand.discard(e)
                cand.add((e[0], e[1], -e[2]))
            else:
                if len(cand) >= max_edges:
                    continue
                pair = full_pairs[rng.integers(len(full_pairs))]
                sign = 1 if rng.random() < 0.5 else -1
                if not _legal_add(cand, pair, sign):
                    continue
                cand.add((pair[0], pair[1], sign))
            if cand != current:
                proposals.append(frozenset(cand))
        if not proposals:
            continue
        scored = score_many(proposals)
        if trace is not None:
            for es, (h, ln) in zip(proposals, scored):
                trace.add(es, h, ln)
        k = int(np.nanargmax([h for h, _ in scored]))
        cand_score = scored[k][0]
        delta = cand_score - cur_score
        if delta > 0 or rng.random() < np.exp(delta / max(temp, 1e-9)):
            current = set(proposals[k])
            cur_score = cand_score
        if cur_score > best_score:
            best, best_score = set(current), cur_score

    return best


def edges_of(spec: ModelSpec) -> set[EdgeKey]:
    return {(e.regulator, e.target, int(e.sign)) for e in spec.edges}
