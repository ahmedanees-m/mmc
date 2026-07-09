"""Logic gates in the grammar: a non-monotone XOR is representable and computes, and
a mutual-repression toggle switch is bistable."""
import numpy as np

from mmc.grammar.model_spec import ModelSpec, Edge, Rule, Term
from mmc.compile.to_ode import build_rhs
from mmc.compile.simulate import steady_state


def _xor():
    spec = ModelSpec(
        genes=["A", "B", "OUT"],
        edges=[Edge(regulator="A", target="OUT", sign=1),
               Edge(regulator="B", target="OUT", sign=1)],
        rules={"OUT": Rule(terms=[
            Term(regulators=["A", "B"], signs={"A": 1, "B": -1}),   # A AND NOT B
            Term(regulators=["A", "B"], signs={"A": -1, "B": 1}),   # NOT A AND B
        ])},
    )
    params = {
        "basal": np.array([0.0, 0.0, 0.0]),
        "decay": np.array([1.0, 1.0, 1.0]),
        "terms": {"OUT": [
            {"prod": 1.0, "w": {"A": 8.0, "B": 8.0}, "theta": {"A": 4.0, "B": -4.0}},
            {"prod": 1.0, "w": {"A": 8.0, "B": 8.0}, "theta": {"A": -4.0, "B": 4.0}},
        ]},
    }
    return spec, params


def test_xor_representable_and_computes():
    spec, params = _xor()
    rhs, idx = build_rhs(spec, params)

    def out(a, b):
        x = np.zeros(3)
        x[idx["A"]], x[idx["B"]] = a, b   # OUT = 0, basal 0, so rhs[OUT] is the production
        return rhs(0.0, x)[idx["OUT"]]

    assert out(1.0, 0.0) > 0.7 and out(0.0, 1.0) > 0.7      # exactly one high -> high
    assert out(0.0, 0.0) < 0.3 and out(1.0, 1.0) < 0.3      # neither or both -> low


def _toggle():
    spec = ModelSpec(
        genes=["A", "B"],
        edges=[Edge(regulator="A", target="B", sign=-1),
               Edge(regulator="B", target="A", sign=-1)],
        rules={"A": Rule(terms=[Term(regulators=["B"])]),
               "B": Rule(terms=[Term(regulators=["A"])])},
    )
    params = {
        "basal": np.array([0.1, 0.1]),
        "decay": np.array([1.0, 1.0]),
        "terms": {
            "A": [{"prod": 4.0, "w": {"B": 8.0}, "theta": -3.0}],
            "B": [{"prod": 4.0, "w": {"A": 8.0}, "theta": -3.0}],
        },
    }
    return spec, params


def test_toggle_switch_is_bistable():
    spec, params = _toggle()
    a_high = steady_state(spec, params, x0=np.array([3.0, 0.0]))
    b_high = steady_state(spec, params, x0=np.array([0.0, 3.0]))
    assert a_high[0] > a_high[1]     # A wins from the A-high basin
    assert b_high[1] > b_high[0]     # B wins from the B-high basin
