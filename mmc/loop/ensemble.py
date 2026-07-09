"""Ensemble tracker.

Steady-state data is degenerate: many structures fit the same training residuals.
Rather than a single winner, keep the structures within a loss margin of the best,
up to a cap, and report the per-edge agreement across them. The frozen artifact is
the ensemble and its edge-agreement profile; the conserved and rewired call in Step 6
is made per edge across the ensemble.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from ..grammar.model_spec import ModelSpec


@dataclass
class Member:
    spec: ModelSpec
    loss: float
    params: dict | None = None


@dataclass
class Ensemble:
    k: int = 5
    margin: float = 0.15
    members: list[Member] = field(default_factory=list)

    def add(self, spec: ModelSpec, loss: float, params: dict | None = None) -> None:
        self.members.append(Member(spec, float(loss), params))
        self.members.sort(key=lambda m: m.loss)
        best = self.members[0].loss
        self.members = [m for m in self.members if m.loss <= best + self.margin][:self.k]

    def best(self) -> Member | None:
        return self.members[0] if self.members else None

    def edge_agreement(self) -> dict[tuple[str, str, int], float]:
        """Fraction of ensemble members that contain each signed edge."""
        n = len(self.members)
        if not n:
            return {}
        counts: Counter = Counter()
        for m in self.members:
            for e in m.spec.edges:
                counts[(e.regulator, e.target, e.sign)] += 1
        return {edge: c / n for edge, c in counts.items()}
