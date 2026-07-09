"""The evaluation prediction and the conserved/rewired decomposition, tested offline."""
import numpy as np

from mmc.eval import conserved_rewired, evaluate
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term


def _spec(edges):
    """Build a spec from (regulator, target, sign) edges with one additive rule each."""
    genes = sorted({g for e in edges for g in (e[0], e[1])})
    regs: dict[str, list[str]] = {}
    for reg, tgt, _ in edges:
        regs.setdefault(tgt, [])
        if reg not in regs[tgt]:
            regs[tgt].append(reg)
    rules = {t: Rule(terms=[Term(regulators=rs)]) for t, rs in regs.items()}
    return ModelSpec(genes=genes,
                     edges=[Edge(regulator=r, target=t, sign=s) for r, t, s in edges],
                     rules=rules)


def test_model_predict_drops_self_and_is_finite():
    spec = _spec([("A", "B", 1), ("B", "C", 1)])
    params = {"basal": np.array([1.0, 0.1, 0.1]), "decay": np.array([1.0, 1.0, 1.0]),
              "terms": {"B": [{"prod": 3.0, "w": {"A": 3.0}, "theta": {"A": 1.0}}],
                        "C": [{"prod": 3.0, "w": {"B": 3.0}, "theta": {"B": 1.0}}]}}
    pred = evaluate.model_predict(spec, params, ["A"], spec.genes)
    assert "A" not in pred["A"]                       # on-target self effect is dropped
    assert all(np.isfinite(v) for v in pred["A"].values())


def test_ge_treats_missing_baseline_as_non_blocking():
    nan = float("nan")
    assert evaluate._ge(0.5, 0.4) and not evaluate._ge(0.4, 0.5)
    assert evaluate._ge(0.5, nan)                     # a degenerate baseline does not block
    assert not evaluate._ge(nan, 0.4)                 # a degenerate model does not pass


def test_decompose_labels_conserved_flipped_dropped_and_added():
    frozen = _spec([("A", "B", 1), ("B", "C", 1), ("A", "C", 1)])
    adapted = _spec([("A", "B", 1), ("B", "C", -1), ("D", "C", 1)])
    dec = conserved_rewired.decompose(frozen, adapted)
    assert ("A", "B") in dec["conserved"]             # same sign across states
    assert ("B", "C") in dec["rewired"]               # sign flip
    assert ("A", "C") in dec["rewired"]               # dropped in adaptation
    assert ("D", "C") in dec["rewired"]               # newly added
