"""Gene-identity anonymisation for the Step 5 ablation (PREREG_v4 section 6).

Step 5 asks a causal question that the literature has so far only answered
correlationally: is the structure the reasoning step proposes conditioned on the
perturbation data, or reconstructed from what it knows about the gene names? The
experiment holds the data fixed and removes the names.

That only works if the names are actually gone. They reach the model through four
separate channels, and missing any one of them defeats the arm silently:

  1. the module gene list in the proposal prompt
  2. the biological context string
  3. the serialised ModelSpec in the repair prompt, which carries gene names in its
     genes list, its edges, and its rule keys
  4. the residual summary, which names perturbations and targets

`redaction_violations` is the audit over the final assembled prompt, and it is what
the test suite asserts on before any arm runs. It is deliberately separate from the
substitution code, so a bug in the substitution cannot also silence its own check.

Substitution is longest-symbol-first with word boundaries. Naive replacement would
turn IL2RB into <alias for IL2>RB, which leaks the identity of both.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from ..grammar.model_spec import Edge, ModelSpec, Rule, Term


@dataclass(frozen=True)
class AliasMap:
    """A bijection between real gene symbols and opaque aliases."""

    to_alias: dict[str, str]
    to_symbol: dict[str, str]

    @property
    def genes(self) -> list[str]:
        return list(self.to_alias)


def make_alias_map(genes: list[str], seed: int = 0, prefix: str = "G") -> AliasMap:
    """Assign aliases in a seed-dependent random order.

    The order is randomised per seed so that no consistent alias survives across runs.
    If GATA3 were always G01, the alias would itself become a name after a few seeds,
    which is the leak this arm exists to close.
    """
    ordered = sorted(set(genes))
    if not ordered:
        return AliasMap({}, {})
    perm = np.random.default_rng(seed).permutation(len(ordered))
    width = max(2, len(str(len(ordered))))
    to_alias = {ordered[int(src)]: f"{prefix}{i + 1:0{width}d}"
                for i, src in enumerate(perm)}
    return AliasMap(to_alias, {a: g for g, a in to_alias.items()})


def _pattern(symbols) -> re.Pattern | None:
    """Longest symbol first, so IL2RB is consumed before IL2 can match inside it."""
    syms = sorted({s for s in symbols if s}, key=len, reverse=True)
    if not syms:
        return None
    return re.compile(r"(?<![A-Za-z0-9_])(" + "|".join(re.escape(s) for s in syms)
                      + r")(?![A-Za-z0-9_])")


def anonymise_text(text: str, amap: AliasMap) -> str:
    """Replace every gene symbol in free text with its alias."""
    pat = _pattern(amap.to_alias)
    if pat is None or not text:
        return text
    return pat.sub(lambda m: amap.to_alias[m.group(1)], text)


def anonymise_spec(spec: ModelSpec, amap: AliasMap) -> ModelSpec:
    """Rewrite a ModelSpec into alias space, structure preserved exactly."""
    a = amap.to_alias
    return ModelSpec(
        genes=[a.get(g, g) for g in spec.genes],
        edges=[Edge(regulator=a.get(e.regulator, e.regulator),
                    target=a.get(e.target, e.target), sign=e.sign) for e in spec.edges],
        rules={a.get(t, t): Rule(terms=[
            Term(regulators=[a.get(r, r) for r in term.regulators],
                 signs={a.get(r, r): s for r, s in term.signs.items()})
            for term in rule.terms]) for t, rule in spec.rules.items()},
    )


def deanonymise_spec(spec: ModelSpec, amap: AliasMap) -> ModelSpec:
    """Map a spec proposed in alias space back to real symbols, for evaluation."""
    return anonymise_spec(spec, AliasMap(amap.to_symbol, amap.to_alias))


def redaction_violations(text: str, genes: list[str]) -> list[str]:
    """Gene symbols still present in an assembled prompt. Empty means the arm is clean.

    Written independently of `anonymise_text` on purpose: if both shared a helper, a
    bug in the matching would hide itself from the audit.
    """
    found = []
    for g in sorted(set(genes)):
        if not g:
            continue
        if re.search(r"(?<![A-Za-z0-9_])" + re.escape(g) + r"(?![A-Za-z0-9_])", text):
            found.append(g)
    return found


def permute_perturbation_labels(observed: np.ndarray, de_mask: np.ndarray,
                                seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Shuffle which response belongs to which perturbation (the A3 and A4 arms).

    Every marginal is preserved and only the pairing between a perturbation and its
    measured response is destroyed, so an arm that scores the same on permuted data
    was not using the pairing. A derangement is used where one exists, since leaving
    a perturbation matched to its own response would weaken the control.
    """
    n = observed.shape[0]
    rng = np.random.default_rng(seed)
    if n < 2:
        return observed.copy(), de_mask.copy()
    for _ in range(100):
        perm = rng.permutation(n)
        if not np.any(perm == np.arange(n)):
            break
    else:
        perm = np.roll(np.arange(n), 1)
    return observed[perm].copy(), de_mask[perm].copy()


def assemble_prompt_surface(genes: list[str], context: str, spec: ModelSpec | None = None,
                            residual_summary: str = "") -> str:
    """Everything the model is shown, concatenated, for the audit to run over.

    Mirrors the four channels in `mmc.shared.llm`. If a channel is added there and not
    here, the redaction test stops covering it, so the two are kept side by side.
    """
    parts = [f"Module genes: {list(genes)}", context or ""]
    if spec is not None:
        parts.append(spec.to_json())
    parts.append(residual_summary or "")
    return "\n\n".join(parts)
