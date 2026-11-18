"""Structural steady-state backend: the fixed point and the interventions."""
import numpy as np

from mmc.compile import structural
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term


def _spec(sign):
    return ModelSpec(
        genes=["A", "B"],
        edges=[Edge(regulator="A", target="B", sign=sign)],
        rules={"B": Rule(terms=[Term(regulators=["A"], signs={"A": sign})])},
    )


def _params():
    return {"basal": [0.5, 0.5],
            "terms": {"B": [{"prod": 2.0, "w": {"A": 3.0}, "theta": {"A": 1.0}}]}}


def test_steady_state_shape_and_finite():
    x = structural.steady_state(_spec(1), _params())
    assert x.shape == (2,)
    assert np.all(np.isfinite(x))


def test_steady_state_respects_clamp():
    x = structural.steady_state(_spec(1), _params(), clamp={"A": 0.0})
    assert x[0] == 0.0


def test_knockdown_of_activator_lowers_target():
    d = structural.knockdown(_spec(1), _params(), "A")
    assert d[0] < 0            # A is driven to zero from its basal level
    assert d[1] < 0            # knocking down an activator lowers its target


def test_activation_of_activator_raises_target():
    d = structural.activation(_spec(1), _params(), "A")
    assert d[1] > 0            # overexpressing an activator raises its target


def test_knockdown_of_repressor_raises_target():
    d = structural.knockdown(_spec(-1), _params(), "A")
    assert d[1] > 0            # knocking down a repressor raises its target


def test_perturb_set_matches_single_knockdown():
    spec, params = _spec(1), _params()
    single = structural.knockdown(spec, params, "A")
    pset = structural.perturb_set(spec, params, ["A"], 0.0)
    assert np.allclose(single, pset)


def test_precomputed_wt_reused():
    spec, params = _spec(1), _params()
    wt = structural.steady_state(spec, params)
    d1 = structural.knockdown(spec, params, "A")
    d2 = structural.knockdown(spec, params, "A", wt=wt)
    assert np.allclose(d1, d2)
