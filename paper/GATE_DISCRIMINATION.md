# Positive control: what the held-out gate does and does not calibrate

*A control for the engineer-behavior calibration claim, measured on the Zhu 2025 CD4+ T-cell
atlas. Scope: this atlas, these edges.*

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
honest and important: the novel hypothesis is **edge-level grounded, not hallucinated** --
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

## Reproduce

`scripts/gate_discrimination.py` builds the hand-specified module, queries the store, and runs
the edge-ablation gate; result in `paper/gate_discrimination.json`.
