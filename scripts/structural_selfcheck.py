"""Backend self-consistency check.

Generate observed knockdown deltas from a known structural model, refit from scratch,
and confirm the JAX structural backend recovers a high in-sample fit. If recovery is
high, the backend is correct and a poor fit on real data is a genuine data or grammar
ceiling, not a bug; if recovery is low, the backend is wrong and must be fixed before any
conclusion.
"""
from __future__ import annotations

import numpy as np

from mmc.compile.structural import knockdown, steady_state
from mmc.eval import metrics
from mmc.fit.fit_structural import multi_fit
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term


def main() -> None:
    genes = ["A", "B", "C", "D", "E"]
    spec = ModelSpec(
        genes=genes,
        edges=[Edge(regulator="A", target="B", sign=1),
               Edge(regulator="B", target="C", sign=-1),
               Edge(regulator="A", target="C", sign=1),
               Edge(regulator="C", target="D", sign=1),
               Edge(regulator="D", target="E", sign=-1)],
        rules={"B": Rule(terms=[Term(regulators=["A"])]),
               "C": Rule(terms=[Term(regulators=["B"]), Term(regulators=["A"])]),
               "D": Rule(terms=[Term(regulators=["C"])]),
               "E": Rule(terms=[Term(regulators=["D"])])},
    )
    true = {"basal": np.array([1.0, 0.2, 0.2, 0.2, 0.2]),
            "terms": {"B": [{"prod": 2.0, "w": {"A": 4.0}, "theta": {"A": 1.0}}],
                      "C": [{"prod": 2.0, "w": {"B": 4.0}, "theta": {"B": 1.0}},
                            {"prod": 1.5, "w": {"A": 3.0}, "theta": {"A": 1.0}}],
                      "D": [{"prod": 2.0, "w": {"C": 4.0}, "theta": {"C": 1.0}}],
                      "E": [{"prod": 2.0, "w": {"D": 4.0}, "theta": {"D": 1.0}}]}}

    gi = {g: i for i, g in enumerate(genes)}
    wt = steady_state(spec, true)
    observed = {reg: {g: float(knockdown(spec, true, reg, wt=wt)[gi[g]])
                      for g in genes if g != reg} for reg in genes}

    fits = multi_fit(spec, observed, n_starts=8, max_iter=400)
    best = fits[0]
    wt2 = steady_state(spec, best["params"])
    pred = {reg: {g: float(knockdown(spec, best["params"], reg, wt=wt2)[gi[g]])
                  for g in genes if g != reg} for reg in genes}
    pooled = metrics.score_set(pred, observed, genes, 0.5)["pooled"]
    print(f"recovery: loss {best['loss']:.4f}  in-sample Pearson {round(pooled['pearson'], 3)}  "
          f"sign-acc {round(pooled['sign_accuracy'], 3)}")


if __name__ == "__main__":
    main()
