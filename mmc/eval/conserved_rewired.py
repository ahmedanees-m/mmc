"""The conserved and rewired decomposition, scored against the Step-0 scaffold.

Compare the structure frozen on the training state with the structure adapted on the
test-state discovery subset. An edge that survives adaptation with the same sign is
called conserved; an edge that flips sign, is dropped, or is added during adaptation is
called rewired. The predicted map is scored against the Stage-0 scaffold, the per-edge
conserved and rewired classification from the precondition test (PREREG section 3),
restricted to the scaffold edges the model has an opinion on. Rewired is the positive
class: precision and recall report how well the loop located the edges that actually
change between the two states.
"""
from __future__ import annotations

from ..data.precondition import classify_module
from ..data.splits import build
from ..grammar.model_spec import ModelSpec


def _edge_signs(spec: ModelSpec) -> dict[tuple[str, str], int]:
    return {(e.regulator, e.target): e.sign for e in spec.edges}


def decompose(frozen: ModelSpec, adapted: ModelSpec) -> dict[str, set]:
    """Predicted conserved and rewired edges from frozen versus adapted structure."""
    f, a = _edge_signs(frozen), _edge_signs(adapted)
    conserved, rewired = set(), set()
    for edge in set(f) | set(a):
        if edge in f and edge in a and f[edge] == a[edge]:
            conserved.add(edge)
        else:
            rewired.add(edge)          # sign flip, dropped, or newly added
    return {"conserved": conserved, "rewired": rewired}


def scaffold(module: str) -> dict[str, set]:
    """Ground-truth conserved and rewired edges from the precondition classification."""
    s = build(module)
    res = classify_module(module, s.train_state, s.test_state)
    conserved = {(e.regulator, e.target) for e in res["edges"] if e.label == "conserved"}
    rewired = {(e.regulator, e.target) for e in res["edges"] if e.label == "rewired"}
    return {"conserved": conserved, "rewired": rewired, "testable": conserved | rewired}


def score(frozen: ModelSpec, adapted: ModelSpec, module: str) -> dict:
    """Score the predicted rewiring map against the scaffold over the edges the model
    covers. Rewired is the positive class."""
    pred = decompose(frozen, adapted)
    truth = scaffold(module)
    model_edges = pred["conserved"] | pred["rewired"]
    evaluated = model_edges & truth["testable"]
    if not evaluated:
        return {"evaluated_edges": 0, "coverage": 0.0,
                "precision": None, "recall": None, "f1": None,
                "predicted_rewired": [], "true_rewired": []}

    pred_rewired = pred["rewired"] & evaluated
    true_rewired = truth["rewired"] & evaluated
    tp = len(pred_rewired & true_rewired)
    fp = len(pred_rewired - true_rewired)
    fn = len(true_rewired - pred_rewired)
    precision = tp / (tp + fp) if (tp + fp) else None
    recall = tp / (tp + fn) if (tp + fn) else None
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall else None)
    return {
        "evaluated_edges": len(evaluated),
        "coverage": len(evaluated) / len(truth["testable"]) if truth["testable"] else 0.0,
        "precision": precision, "recall": recall, "f1": f1,
        "n_predicted_rewired": len(pred_rewired), "n_true_rewired": len(true_rewired),
        "predicted_rewired": sorted(pred_rewired), "true_rewired": sorted(true_rewired),
    }
