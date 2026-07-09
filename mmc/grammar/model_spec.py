"""The executable-model grammar.

The reasoning step emits structure: genes, signed edges, and rules (the logic
form). The optimizer (fit/) sets magnitudes (params). This split is the actor and
tool boundary.

Rules are a bounded sum-of-products over sigmoid gates (interpretable DNF):
    production_i = sum_k prod_ik * PROD_{j in term_k} sigma(sign*w_ikj * x_j - theta_ik)
A single additive term = the monotone default. Product of sigmoids = AND.
Sum of terms = OR. Negative weight = NOT. OR-of-ANDs = XOR / non-monotone logic.
Bounded (<=3 terms/target, <=3 regulators/term) to preserve interpretability
and identifiability -- NOT arbitrary MLPs.
"""
from __future__ import annotations
from pydantic import BaseModel, field_validator, model_validator

MAX_TERMS = 3
MAX_REGS_PER_TERM = 3


class Edge(BaseModel):
    regulator: str
    target: str
    sign: int  # +1 activation, -1 repression

    @field_validator("sign")
    @classmethod
    def _sign(cls, v: int) -> int:
        if v not in (1, -1):
            raise ValueError("sign must be +1 or -1")
        return v


class Term(BaseModel):
    """One product-of-sigmoids gate over a small set of regulators."""
    regulators: list[str]

    @field_validator("regulators")
    @classmethod
    def _size(cls, v: list[str]) -> list[str]:
        if not (1 <= len(v) <= MAX_REGS_PER_TERM):
            raise ValueError(f"term needs 1..{MAX_REGS_PER_TERM} regulators")
        return v


class Rule(BaseModel):
    """Sum of product-terms for one target (DNF over gates)."""
    terms: list[Term]

    @field_validator("terms")
    @classmethod
    def _nterms(cls, v: list[Term]) -> list[Term]:
        if not (1 <= len(v) <= MAX_TERMS):
            raise ValueError(f"rule needs 1..{MAX_TERMS} terms")
        return v

    @property
    def is_additive(self) -> bool:
        return len(self.terms) == 1 and len(self.terms[0].regulators) >= 1


class ModelSpec(BaseModel):
    genes: list[str]
    edges: list[Edge]
    rules: dict[str, Rule]  # keyed by target gene

    @model_validator(mode="after")
    def _coherent(self) -> "ModelSpec":
        gs = set(self.genes)
        if len(gs) != len(self.genes):
            raise ValueError("duplicate genes")
        for e in self.edges:
            if e.regulator not in gs or e.target not in gs:
                raise ValueError(f"dangling edge {e.regulator}->{e.target}")
        in_edges: dict[str, set[str]] = {g: set() for g in self.genes}
        for e in self.edges:
            in_edges[e.target].add(e.regulator)
        for tgt, rule in self.rules.items():
            if tgt not in gs:
                raise ValueError(f"rule for unknown gene {tgt}")
            for term in rule.terms:
                for r in term.regulators:
                    if r not in in_edges[tgt]:
                        raise ValueError(
                            f"rule term for {tgt} uses {r} without an edge {r}->{tgt}"
                        )
        return self

    # round-trip helpers
    def to_json(self) -> str:
        return self.model_dump_json()

    @classmethod
    def from_json(cls, s: str) -> "ModelSpec":
        return cls.model_validate_json(s)

    def edge_sign(self, regulator: str, target: str) -> int | None:
        for e in self.edges:
            if e.regulator == regulator and e.target == target:
                return e.sign
        return None
