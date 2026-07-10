"""Simulator core: WT steady state and knockdown direction."""
import numpy as np

from mmc.compile.perturb import knockdown
from mmc.compile.simulate import steady_state
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term


def _model():
    spec = ModelSpec(
        genes=["A", "B"],
        edges=[Edge(regulator="A", target="B", sign=1)],
        rules={"B": Rule(terms=[Term(regulators=["A"])])},
    )
    params = {
        "basal": np.array([1.0, 0.05]),   # A is a source at ~1.0
        "decay": np.array([1.0, 1.0]),
        "terms": {"B": [{"prod": 3.0, "w": {"A": 3.0}, "theta": 1.0}]},
    }
    return spec, params


def test_wt_positive():
    spec, params = _model()
    ss = steady_state(spec, params)
    assert ss.shape == (2,)
    assert ss[0] > 0.5           # A source
    assert ss[1] > ss[0]         # B activated above A's level


def test_knockdown_lowers_activator_and_target():
    spec, params = _model()
    d = knockdown(spec, params, "A")
    assert d[0] < 0              # A knocked down
    assert d[1] < 0              # B loses its activator
