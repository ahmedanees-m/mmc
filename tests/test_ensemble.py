"""The ensemble tracker: it keeps structures within a loss margin, up to a cap, and
reports per-edge agreement."""
from mmc.grammar.model_spec import ModelSpec, Edge, Rule, Term
from mmc.loop.ensemble import Ensemble


def _spec(sign: int):
    return ModelSpec(
        genes=["A", "B"], edges=[Edge(regulator="A", target="B", sign=sign)],
        rules={"B": Rule(terms=[Term(regulators=["A"])])},
    )


def test_keeps_within_margin_and_reports_best():
    e = Ensemble(k=5, margin=0.1)
    e.add(_spec(1), 1.00)
    e.add(_spec(-1), 1.05)
    e.add(_spec(1), 2.00)          # outside the margin, dropped
    assert len(e.members) == 2
    assert e.best().loss == 1.00


def test_edge_agreement_fractions():
    e = Ensemble(k=5, margin=1.0)
    e.add(_spec(1), 1.0)
    e.add(_spec(1), 1.1)
    e.add(_spec(-1), 1.2)
    agreement = e.edge_agreement()
    assert abs(agreement[("A", "B", 1)] - 2 / 3) < 1e-9
    assert abs(agreement[("A", "B", -1)] - 1 / 3) < 1e-9
