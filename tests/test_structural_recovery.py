"""The structural backend recovers a known generating model.

Observed knockdown deltas are generated from a model whose parameters are known, the
parameters are then refit from scratch, and the refitted predictions are compared against
the generated ones. A low correlation here would mean a poor fit on real data could be a
backend fault rather than a property of the data, so this is the check that separates the
two.
"""
import numpy as np
import pytest

pytest.importorskip("jax", reason="the structural fit path needs the JAX backend")

from mmc.compile.structural import knockdown, steady_state  # noqa: E402
from mmc.fit.fit_structural import multi_fit  # noqa: E402
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term  # noqa: E402

GENES = ["A", "B", "C", "D", "E"]

SPEC = ModelSpec(
    genes=GENES,
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

TRUE = {"basal": np.array([1.0, 0.2, 0.2, 0.2, 0.2]),
        "terms": {"B": [{"prod": 2.0, "w": {"A": 4.0}, "theta": {"A": 1.0}}],
                  "C": [{"prod": 2.0, "w": {"B": 4.0}, "theta": {"B": 1.0}},
                        {"prod": 1.5, "w": {"A": 3.0}, "theta": {"A": 1.0}}],
                  "D": [{"prod": 2.0, "w": {"C": 4.0}, "theta": {"C": 1.0}}],
                  "E": [{"prod": 2.0, "w": {"D": 4.0}, "theta": {"D": 1.0}}]}}


def _responses(params):
    """Every single-gene knockdown's effect on every other gene, as a flat vector."""
    gi = {g: i for i, g in enumerate(GENES)}
    wt = steady_state(SPEC, params)
    return {reg: {g: float(knockdown(SPEC, params, reg, wt=wt)[gi[g]])
                  for g in GENES if g != reg} for reg in GENES}


def _flatten(d):
    return np.array([d[reg][g] for reg in GENES for g in GENES if g != reg], float)


def test_the_backend_recovers_a_known_generating_model():
    observed = _responses(TRUE)
    best = multi_fit(SPEC, observed, n_starts=8, max_iter=400)[0]
    pred = _responses(best["params"])

    o, p = _flatten(observed), _flatten(pred)
    pearson = float(np.corrcoef(o, p)[0, 1])

    assert pearson == pytest.approx(1.0, abs=1e-4), f"in-sample Pearson {pearson:.6f}"
    assert np.max(np.abs(o - p)) < 1e-3

    # signs agree on every response large enough for a sign to mean anything; the
    # remainder are zeros, where the fitted value differs only in the sign of the zero
    big = np.abs(o) > 1e-6
    assert np.all(np.sign(o[big]) == np.sign(p[big]))
