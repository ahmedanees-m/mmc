"""Verify the CRISPRa activation operator on a toy (PREREG_norman.md section 3).

The Norman perturbation is do(x = high), not the do(x = 0) knockdown used for Zhu. On a
known circuit, activating an activator must raise its targets, the opposite of knocking it
down, and a double must clamp both genes. If these hold, the operator is correct and the
Norman fit can proceed; if not, the whole Norman run would be invalid.
"""
from __future__ import annotations

import numpy as np

from mmc.compile.structural import activation, knockdown, perturb_set, steady_state
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term

spec = ModelSpec(
    genes=["A", "B", "C", "D"],
    edges=[Edge(regulator="A", target="B", sign=1),
           Edge(regulator="B", target="C", sign=1),
           Edge(regulator="A", target="D", sign=-1)],
    rules={"B": Rule(terms=[Term(regulators=["A"])]),
           "C": Rule(terms=[Term(regulators=["B"])]),
           "D": Rule(terms=[Term(regulators=["A"])])},
)
params = {"basal": np.array([1.0, 0.1, 0.1, 1.0]),
          "terms": {"B": [{"prod": 3.0, "w": {"A": 3.0}, "theta": {"A": 1.0}}],
                    "C": [{"prod": 3.0, "w": {"B": 3.0}, "theta": {"B": 1.0}}],
                    "D": [{"prod": 3.0, "w": {"A": 3.0}, "theta": {"A": 1.0}}]}}

gi = {g: i for i, g in enumerate(spec.genes)}
wt = steady_state(spec, params)

act = activation(spec, params, "A", level=4.0, wt=wt)     # A activates B and C, represses D
kd = knockdown(spec, params, "A", wt=wt)
print(f"activate A: dB {act[gi['B']]:+.2f}  dC {act[gi['C']]:+.2f}  dD {act[gi['D']]:+.2f}")
print(f"knockdown A: dB {kd[gi['B']]:+.2f}  dC {kd[gi['C']]:+.2f}  dD {kd[gi['D']]:+.2f}")

assert act[gi["B"]] > 0 and act[gi["C"]] > 0, "activating an activator must raise its targets"
assert act[gi["D"]] < 0, "activating a repressor's source must lower the target"
assert kd[gi["B"]] < 0 and kd[gi["C"]] < 0, "knockdown must lower an activator's targets"
assert (np.sign(act[gi["B"]]) != np.sign(kd[gi["B"]])), "activation and knockdown must oppose"

double = perturb_set(spec, params, ["A", "B"], level=4.0, wt=wt)   # clamp both high
assert double[gi["C"]] > 0, "double activation propagates to C"
print(f"double activate A,B: dC {double[gi['C']]:+.2f}")
print("[toycheck] activation operator do(x=high) is correct.")
