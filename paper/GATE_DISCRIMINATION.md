# Positive control: what the held-out gate does and does not calibrate

*A control for the engineer-behavior calibration claim, measured on the Zhu 2025 CD4+ T-cell
atlas. Scope: this atlas, these edges.*

**The finding.** The failure mode of AI-driven mechanistic discovery here is not hallucination.
The loop's novel hypotheses are individually grounded and pass edge-level necessity and ablation
exactly like textbook biology; what they fail is composing into a model that beats a strong
baseline held-out. Grounded-but-non-predictive hypotheses pass the interpretability checks
(necessity, ablation) that most validation relies on, and the only safeguard that catches them
is held-out predictive advantage over a strong baseline. This control establishes both halves:
the module-level gate discriminates (it says yes to a true structure and ranks an over-connected
one below it), and the novel edges pass edge-level necessity yet the model still does not beat
the baseline.

## Why this control

The engineer-behavior result is that the loop's plausible novel hypotheses do not make the
model beat a linear baseline held-out. On its own that does not show the held-out gate is a
*discriminator* rather than a blanket rejecter: a gate that says no to everything is not
demonstrated to say yes to anything true. So the control asks two questions of the per-edge
gate (an edge is REQUIRED only if removing it lowers the model's own held-out ACC_DEG):

1. Does it ever say yes to a true edge? (the positive control)
2. Does it discriminate a true edge from a plausible-but-novel one?

## Method

One strong-signal module holds the textbook Th2 circuit (GATA3 -> IL5 is the single strongest
knockdown effect in the atlas) together with the demo's plausible-but-novel STK11 -> chemokine
hypothesis. All 11 genes are perturbed in the atlas, giving 11 leave-one-out folds. The
structure is hand-built (no model calls); only the store is queried. The same edge-ablation
gate is run on every edge (`scripts/gate_discrimination.py`).

## Result

Held-out ACC_DEG drop when the edge is removed (base held-out ACC_DEG 0.5):

| Edge | Class | Held-out drop | Required |
|---|---|---|---|
| GATA3 -> IL13 | textbook | 0.292 | yes |
| STK11 -> CCL3 | novel | 0.250 | yes |
| STK11 -> CCL4 | novel | 0.250 | yes |
| STAT6 -> GATA3 | textbook | 0.208 | yes |
| GATA3 -> IL5 | textbook | 0.167 | yes |
| STK11 -> CXCL8 | novel | 0.125 | yes |
| GATA3 -> IL4 | textbook | 0.083 | yes |
| STAT4 -> TBX21 | textbook | 0.083 | yes |
| TBX21 -> GATA3 | textbook | 0.083 | yes |

## What it means

**The positive control passes.** The gate says yes to true edges: removing the textbook
GATA3 -> IL5, GATA3 -> IL13, or STAT6 -> GATA3 edges lowers the model's held-out accuracy. The
gate is not a blanket rejecter.

**But the edge-ablation gate does not discriminate textbook from novel.** The STK11 -> chemokine
edges are flagged required too, some with a larger drop than the textbook edges. This is
important: the novel hypothesis is **edge-level grounded, not a spurious edge**, because
STK11 knockdown genuinely moves those chemokines. Edge-level necessity is a weak bar: it
largely confirms that the edge's regulator, when knocked down, moves its target, which is what
the edge was drawn from. It does not test advantage over a baseline.

**The calibration that matters is the module-level held-out gate, not edge ablation.** The bar
that separates a real marginal effect from a predictive model is whether the mechanistic model
beats a linear baseline on held-out prediction. Every module fails it (0 of 2; the cytokine
module 0.18 versus linear 0.45, separated CIs). So the individually-grounded edges do not
compose into a model that beats a simple baseline, and it is the module-level gate that refuses
to certify the STK11-containing model.

The corrected lesson, scoped precisely: **plausible mechanism, even when each edge is grounded
in a real marginal effect, does not equal a model that predicts better than a simple baseline.**
The engineer's hypotheses are not hallucinations; they are real effects that do not compose
into predictive advantage, and the module-level held-out gate is what reveals the gap. This is
a more precise claim than "the AI's plausible ideas are wrong," and it is the one the data
supports.

## The module-level gate is a verified discriminator, not a uniform rejecter

The corrected claim rests on the module-level held-out-advantage test, so a synthetic control
verifies that the test can fire positive rather than rejecting everything. On synthetic ground truth (data generated
from a known sparse structure with additive noise, `scripts/gate_synthetic_control.py`), the
true structure clears the held-out bar cleanly: leave-one-perturbation-out DE-overlap 0.90
(95% CI [0.70, 1.00]) and ACC_DEG 1.0, against the mean baseline's 0.17 and 0.37. A
fully-connected over-connected structure over-fits and drops to DE-overlap 0.63 and ACC_DEG
0.77, well below the true model. So the test certifies a true structure and ranks an
over-connected one materially lower. On every real module it rejects, including the grounded
STK11 model, it does so because none beats a strong linear baseline, not because it rejects
everything. Edge-level necessity is the weaker bar that the grounded edges pass.

## Reproduce

`scripts/gate_discrimination.py` builds the hand-specified module, queries the store, and runs
the edge-ablation gate (`paper/gate_discrimination.json`). `scripts/gate_synthetic_control.py`
runs the synthetic discriminator control (`paper/gate_synthetic_control.json`).
