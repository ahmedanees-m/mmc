"""Residual reader: a compact, structured summary of the systematic structural
failures, for the repair step.

Only structural residuals (from the diagnostic gate) are read; parametric ones are
withheld, since they are the optimizer's to resolve. Each is grouped by pattern:
wrong-sign (the model moves the target the wrong way), missing-effect (the model
predicts no effect where the data has a large one), or spurious-effect (the model
predicts a large effect where the data has none). The summary is bounded to the
largest gaps so it stays within a token budget.
"""
from __future__ import annotations

from ..fit import diagnose as diag


def _pattern(pred: float, obs: float) -> str:
    if abs(obs) > 0.5 and abs(pred) > 0.5 and (pred > 0) != (obs > 0):
        return "wrong_sign"
    if abs(pred) < 0.3 and abs(obs) > 0.7:
        return "missing_effect"
    if abs(pred) > 0.7 and abs(obs) < 0.3:
        return "spurious_effect"
    return "magnitude"


def structural_items(spec, best_fit: dict, labels: dict, observed: dict,
                     top: int = 12, residuals_fn=None) -> list[dict]:
    """The largest structural residuals, each with its failure pattern. residuals_fn
    selects the backend that recomputes predictions (defaults to the ODE backend)."""
    res = (residuals_fn or diag.residuals)(spec, best_fit["params"], observed)
    items = []
    for key, (pred, obs) in res.items():
        if labels.get(key) != "structural":
            continue
        pert, target = key
        items.append({"perturbation": pert, "target": target,
                      "predicted": round(float(pred), 2), "observed": round(float(obs), 2),
                      "gap": round(abs(float(obs) - float(pred)), 2),
                      "pattern": _pattern(float(pred), float(obs))})
    items.sort(key=lambda d: -d["gap"])
    return items[:top]


def summary_text(items: list[dict], stats: dict | None = None) -> str:
    if not items:
        return "No structural residuals remain; the current structure is consistent with the data."
    lines = ["Structural residuals (knockdown, target: model delta vs observed delta):"]
    for d in items:
        lines.append(
            f"- {d['perturbation']} knockdown on {d['target']}: model {d['predicted']:+.2f}, "
            f"observed {d['observed']:+.2f} ({d['pattern']})")
    if stats is not None:
        lines.append(f"Fit confidence: best loss {stats.get('loss_best', float('nan')):.3f}, "
                     f"seed spread {stats.get('loss_spread', float('nan')):.3f} over "
                     f"{stats.get('n_starts', 0)} starts.")
    return "\n".join(lines)
