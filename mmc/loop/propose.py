"""Propose an initial structure via the model client, logging to the trace.

The reasoning step sees only the module gene list, the biological context, the
current ModelSpec, and the structural residual summary. It never sees held-out or
test-state data (Tier A) or the Tier-B held-out subset. Every proposal and its
rationale is logged to the trace.
"""
from __future__ import annotations

import hashlib

from ..grammar.model_spec import ModelSpec
from ..shared import llm
from ..shared.trace import Step, Trace


def model_hash(spec: ModelSpec) -> str:
    return hashlib.sha256(spec.to_json().encode()).hexdigest()[:10]


def propose(module_genes: list[str], context: str,
            trace: Trace | None = None) -> tuple[ModelSpec, str]:
    spec, rationale = llm.propose_structure(module_genes, context)
    if trace is not None:
        trace.log(Step(kind="propose", model_hash=model_hash(spec),
                       rationale=rationale or "initial structure",
                       edit={"edges": len(spec.edges)}, train_score=None))
    return spec, rationale
