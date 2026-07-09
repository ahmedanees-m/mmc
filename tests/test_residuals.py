"""The residual reader: pattern classification, summary text, and structural
filtering, tested offline (no model client, no cma)."""
import numpy as np

from mmc.grammar.model_spec import ModelSpec, Edge, Rule, Term
from mmc.loop.residuals import _pattern, structural_items, summary_text


def test_pattern_classification():
    assert _pattern(-2.0, +2.0) == "wrong_sign"
    assert _pattern(0.0, +2.0) == "missing_effect"
    assert _pattern(+2.0, 0.0) == "spurious_effect"


def test_summary_text_reports_no_residuals():
    assert "No structural residuals" in summary_text([])


def test_summary_text_lists_items():
    items = [{"perturbation": "A", "target": "B", "predicted": -1.0,
              "observed": 2.0, "gap": 3.0, "pattern": "wrong_sign"}]
    text = summary_text(items)
    assert "A knockdown on B" in text and "wrong_sign" in text


def test_structural_items_filters_to_structural_and_labels_pattern():
    spec = ModelSpec(
        genes=["A", "B"], edges=[Edge(regulator="A", target="B", sign=1)],
        rules={"B": Rule(terms=[Term(regulators=["A"])])},
    )
    params = {"basal": np.array([1.0, 0.05]), "decay": np.array([1.0, 1.0]),
              "terms": {"B": [{"prod": 3.0, "w": {"A": 3.0}, "theta": 1.0}]}}
    best = {"params": params}
    observed = {"A": {"B": +2.0}}                 # B rises on A knockdown: wrong sign
    labels = {("A", "B"): "structural"}
    items = structural_items(spec, best, labels, observed)
    assert len(items) == 1
    assert items[0]["pattern"] == "wrong_sign"
