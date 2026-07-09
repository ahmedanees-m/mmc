"""Grammar schema round-trip and validators."""
import pytest
from mmc.grammar.model_spec import ModelSpec, Edge, Rule, Term


def _spec():
    return ModelSpec(
        genes=["A", "B"],
        edges=[Edge(regulator="A", target="B", sign=1)],
        rules={"B": Rule(terms=[Term(regulators=["A"])])},
    )


def test_round_trip():
    s = _spec()
    s2 = ModelSpec.from_json(s.to_json())
    assert s2.genes == s.genes
    assert s2.edge_sign("A", "B") == 1


def test_dangling_edge_rejected():
    with pytest.raises(Exception):
        ModelSpec(genes=["A"], edges=[Edge(regulator="A", target="Z", sign=1)], rules={})


def test_rule_without_edge_rejected():
    with pytest.raises(Exception):
        ModelSpec(
            genes=["A", "B"], edges=[],
            rules={"B": Rule(terms=[Term(regulators=["A"])])},
        )


def test_bad_sign_rejected():
    with pytest.raises(Exception):
        Edge(regulator="A", target="B", sign=0)


def test_signed_regulator_form_is_normalized():
    # a term's regulators given as objects with signs is split into names plus a signs map
    spec = ModelSpec(
        genes=["A", "B"],
        edges=[Edge(regulator="A", target="B", sign=1)],
        rules={"B": {"terms": [{"regulators": [{"regulator": "A", "sign": -1}]}]}},
    )
    term = spec.rules["B"].terms[0]
    assert term.regulators == ["A"]
    assert term.signs == {"A": -1}
