"""The model-evolution trace: every proposal, edit, rationale, and score.

The trace is the auditable record of how the model reasoned across the loop,
including hypotheses that were later falsified against data. The pattern is
carried over from VERDICT's inspectable trace.
"""
from __future__ import annotations
import json, time
from dataclasses import dataclass, field, asdict
from pathlib import Path


@dataclass
class Step:
    kind: str                 # 'propose' | 'repair' | 'freeze'
    model_hash: str
    rationale: str            # the reasoning step's mechanistic justification
    edit: dict | None         # the structural change, if any
    train_score: float | None
    ts: float = field(default_factory=time.time)


@dataclass
class Trace:
    module: str
    steps: list[Step] = field(default_factory=list)

    def log(self, step: Step) -> None:
        self.steps.append(step)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps({"module": self.module,
                        "steps": [asdict(s) for s in self.steps]}, indent=2)
        )
