"""v3 Part 2d: the known-structure probe.

Fit the textbook CD4 T-helper transcription-factor network under the structural
steady-state backend and report the in-sample fit and the wall time. A random-structure
control fits the same way. This resolves the Part-0 gate: if the textbook structure fits
materially better than the ODE ceiling (Pearson about 0.30) and better than random, the
earlier failure was a search and backend problem, the engine works, and discovery can
proceed; if it fits no better, the ceiling is the grammar or the data, and discovery plus
the limit-map carry the paper.
"""
from __future__ import annotations

import time

import numpy as np

from mmc.compile.structural import knockdown, steady_state
from mmc.data import module_extract
from mmc.eval import metrics
from mmc.fit.fit_structural import multi_fit
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term

# Well-established CD4 lineage TF regulatory edges (regulator, target, sign). Additive:
# each regulator is its own single-regulator term, so a target's regulators contribute
# independently, up to the grammar's three-term bound.
TEXTBOOK = [
    ("STAT6", "GATA3", 1), ("TBX21", "GATA3", -1), ("RUNX3", "GATA3", -1),
    ("STAT4", "TBX21", 1), ("STAT1", "TBX21", 1), ("GATA3", "TBX21", -1),
    ("STAT3", "RORC", 1), ("FOXP3", "RORC", -1), ("BATF", "RORC", 1),
    ("RORC", "FOXP3", -1),
    ("STAT3", "BCL6", 1), ("PRDM1", "BCL6", -1),
    ("BCL6", "PRDM1", -1), ("IRF4", "PRDM1", 1),
    ("IRF4", "BATF", 1),
    ("IRF4", "MAF", 1), ("BATF", "MAF", 1),
]


def _spec_from_edges(genes: list[str], triples: list[tuple[str, str, int]]) -> ModelSpec:
    triples = [(r, t, s) for r, t, s in triples if r in genes and t in genes and r != t]
    regs_by_tgt: dict[str, list[tuple[str, int]]] = {}
    for r, t, s in triples:
        regs_by_tgt.setdefault(t, [])
        if r not in [x for x, _ in regs_by_tgt[t]]:
            regs_by_tgt[t].append((r, s))
    edges, rules = [], {}
    for t, rs in regs_by_tgt.items():
        rs = rs[:3]                                    # at most three terms per target
        for r, s in rs:
            edges.append(Edge(regulator=r, target=t, sign=s))
        rules[t] = Rule(terms=[Term(regulators=[r]) for r, _ in rs])
    return ModelSpec(genes=list(genes), edges=edges, rules=rules)


def _random_spec(genes: list[str], n_edges: int, seed: int) -> ModelSpec:
    rng = np.random.default_rng(seed)
    triples, seen = [], set()
    while len(triples) < n_edges:
        r, t = rng.choice(genes, 2, replace=False)
        if (r, t) in seen:
            continue
        seen.add((r, t))
        triples.append((str(r), str(t), int(rng.choice([1, -1]))))
    return _spec_from_edges(genes, triples)


def _predict(spec, params, perts, genes):
    gi = {g: i for i, g in enumerate(spec.genes)}
    wt = steady_state(spec, params)
    out = {}
    for p in perts:
        if p in gi:
            d = knockdown(spec, params, p, wt=wt)
            out[p] = {g: float(d[gi[g]]) for g in genes if g in gi and g != p}
    return out


def _evaluate(name, spec, observed, genes, n_starts, max_iter):
    t0 = time.time()
    fits = multi_fit(spec, observed, n_starts=n_starts, max_iter=max_iter)
    dt = time.time() - t0
    best = fits[0]
    pred = _predict(spec, best["params"], list(observed), genes)
    pooled = metrics.score_set(pred, observed, genes, 0.5)["pooled"]
    r, s = pooled["pearson"], pooled["sign_accuracy"]
    print(f"  {name:22s} edges {len(spec.edges):>3d}  loss {best['loss']:.4f}  "
          f"in-sample Pearson {r if r is None else round(r, 3)}  "
          f"sign-acc {s if s is None else round(s, 3)}  [{dt:.1f}s, {n_starts} starts]")


def main() -> None:
    genes = module_extract.model_genes("CD4_TF_network")
    observed = module_extract.observed_deltas("CD4_TF_network", "Stim8hr")
    print(f"known-structure probe on CD4_TF_network: {len(genes)} genes, "
          f"{len(observed)} perturbations, structural backend")
    textbook = _spec_from_edges(genes, TEXTBOOK)
    _evaluate("textbook", textbook, observed, genes, n_starts=10, max_iter=300)
    for seed in (0, 1, 2):
        _evaluate(f"random (seed {seed})", _random_spec(genes, len(textbook.edges), seed),
                  observed, genes, n_starts=10, max_iter=300)


if __name__ == "__main__":
    main()
