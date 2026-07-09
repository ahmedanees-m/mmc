"""Reasoning-guided structural repair via the model client, logging to the trace.

The reasoning step reads the structured residual summary and infers a mechanism: a
knockdown that raises a target the model lowers implies a sign flip or an
intermediate repressor; a residual that is non-monotone in two regulators implies
replacing an additive term with a product-term gate. It emits a structural edit. The
loop re-fits, re-diagnoses, re-reads, and iterates to a plateau or a budget.
"""
from __future__ import annotations

from ..grammar.model_spec import ModelSpec
from ..shared import llm
from ..shared.trace import Step, Trace
from .propose import model_hash


def _diff(old: ModelSpec, new: ModelSpec) -> dict:
    old_e = {(e.regulator, e.target, e.sign) for e in old.edges}
    new_e = {(e.regulator, e.target, e.sign) for e in new.edges}
    return {"added": sorted(f"{r}->{t}({s:+d})" for r, t, s in new_e - old_e),
            "removed": sorted(f"{r}->{t}({s:+d})" for r, t, s in old_e - new_e)}


def repair(current: ModelSpec, context: str, residual_summary: str,
           trace: Trace | None = None) -> tuple[ModelSpec, str]:
    spec, rationale = llm.repair_structure(current, context, residual_summary)
    if trace is not None:
        trace.log(Step(kind="repair", model_hash=model_hash(spec),
                       rationale=rationale or "structural edit",
                       edit=_diff(current, spec), train_score=None))
    return spec, rationale
