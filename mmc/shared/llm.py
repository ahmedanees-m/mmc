"""Model client for the reasoning step.

Carried over from the VERDICT engine client: structured output, adaptive thinking,
an effort control, and refusal handling. In MMC the reasoning step proposes and
repairs a ModelSpec (grammar.model_spec); it sees only the module gene list, the
biological context, the current structure, and the structural residual summary. It
never sees held-out or test-state data. The reasoning step sets structure and logic
form; the optimizer sets every magnitude.

The key is read from the environment; it is never written to the repository. The
model is configurable with MMC_MODEL and defaults to claude-opus-4-8.
"""
from __future__ import annotations

import json
import os
import re
from functools import lru_cache

from ..grammar.model_spec import ModelSpec

DEFAULT_MODEL = "claude-opus-4-8"


def model() -> str:
    return os.environ.get("MMC_MODEL", DEFAULT_MODEL)


def _supports_effort(m: str) -> bool:
    """Adaptive thinking and the effort control are available on Opus 4.6 and later,
    Sonnet 4.6 and later, and Fable 5. Haiku and older tiers take neither."""
    m = m.lower()
    if "haiku" in m:
        return False
    return m not in {"claude-sonnet-4-5", "claude-opus-4-1", "claude-opus-4-0"}


@lru_cache(maxsize=1)
def _client():
    import anthropic  # imported lazily so the package loads without the SDK configured
    return anthropic.Anthropic()


def _text(resp) -> str:
    if getattr(resp, "stop_reason", None) == "refusal":
        raise RuntimeError("The model declined the request (stop_reason=refusal).")
    return "".join(b.text for b in resp.content if b.type == "text").strip()


def reason(system: str, user: str, *, effort: str = "high", max_tokens: int = 8000) -> str:
    """A reasoning call: adaptive thinking on where supported, plain text out."""
    m = model()
    kwargs = dict(
        model=m, max_tokens=max_tokens, system=system,
        messages=[{"role": "user", "content": user}],
    )
    if _supports_effort(m):
        kwargs["thinking"] = {"type": "adaptive"}
        kwargs["output_config"] = {"effort": effort}
    return _text(_client().messages.create(**kwargs))


def extract(system: str, user: str, schema: dict, *, max_tokens: int = 2000) -> dict:
    """An extraction call: output constrained to `schema`, no thinking."""
    resp = _client().messages.create(
        model=model(), max_tokens=max_tokens, system=system,
        output_config={"format": {"type": "json_schema", "schema": schema}},
        messages=[{"role": "user", "content": user}],
    )
    return json.loads(_text(resp))


# ---- structure proposal and repair over the grammar ----

_SYSTEM = (
    "You are the structure proposer for a mechanistic model-discovery loop over CD4+ "
    "T-cell gene regulation. Given a module gene list and biological context, express a "
    "gene-regulatory circuit as a ModelSpec: genes, signed edges (regulator, target, sign "
    "+1 or -1), and rules (a bounded sum of product-terms over sigmoid gates for each "
    "target). Default each rule to a single additive term; use a product-term gate only "
    "when non-monotone logic is required. You set structure and logic form only; an "
    "optimizer sets all magnitudes, so do not emit parameter values. Every regulator named "
    "in a rule term must have a corresponding edge into that target. A term's regulators "
    "is a list of gene-name strings; put any per-term gate sign in an optional signs "
    "object mapping a regulator to +1 or -1. Example: "
    '{"genes": ["A", "B", "C"], '
    '"edges": [{"regulator": "A", "target": "C", "sign": 1}, '
    '{"regulator": "B", "target": "C", "sign": -1}], '
    '"rules": {"C": {"terms": [{"regulators": ["A", "B"], "signs": {"B": -1}}]}}}. '
    "Begin with a single line starting 'Rationale:' giving one sentence on the key "
    "structural choice, then emit the ModelSpec as one JSON object inside a single fenced "
    "json block and nothing else."
)


def _json_block(text: str) -> str:
    m = re.search(r"```(?:json)?\s*(\{.*\})\s*```", text, re.S) or re.search(r"(\{.*\})", text, re.S)
    if not m:
        raise ValueError("no JSON object found in the model output")
    return m.group(1)


def _rationale(text: str) -> str:
    m = re.search(r"Rationale:\s*(.+)", text)
    return m.group(1).strip() if m else ""


def _emit(user: str, *, attempts: int = 2) -> tuple[ModelSpec, str]:
    prompt = user
    last = ""
    for _ in range(attempts):
        out = reason(_SYSTEM, prompt, effort="high")
        try:
            return ModelSpec.from_json(_json_block(out)), _rationale(out)
        except Exception as e:  # invalid schema or JSON; return the reason and retry
            last = str(e)
            prompt = user + (f"\n\nYour previous output was rejected: {last}. "
                             "Re-emit one valid ModelSpec JSON object only.")
    raise RuntimeError(f"structure not valid after {attempts} attempts: {last}")


def propose_structure(module_genes: list[str], context: str) -> tuple[ModelSpec, str]:
    """Propose an initial executable model. Returns (spec, rationale)."""
    user = f"Module genes: {module_genes}\n\nBiological context:\n{context}\n"
    return _emit(user)


def repair_structure(current: ModelSpec, context: str,
                     residual_summary: str) -> tuple[ModelSpec, str]:
    """Revise the structure to address the structural residuals. Returns (spec, rationale)."""
    user = (
        f"Module genes: {current.genes}\n\nBiological context:\n{context}\n\n"
        f"Current model (ModelSpec JSON):\n{current.to_json()}\n\n"
        f"Structural residuals to address:\n{residual_summary}\n\n"
        "Propose the smallest structural edit that addresses them (flip a sign, add or "
        "remove an edge, add a regulator, or escalate a rule to a product-term gate)."
    )
    return _emit(user)
