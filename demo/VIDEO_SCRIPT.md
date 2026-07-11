# Demo walkthrough (3 minutes)

A narration and screen-capture guide for a short walkthrough of `demo/app.py`. Record the
narration and the screen capture separately and cut to 3:00 or under. Each section maps to one
view in the app.

---

**1. Interrogate the circuit (0:00-1:05)**

*On screen:* the "Interrogate the circuit" view. Select GATA3 knockdown, then STAT6 knockdown,
then the GATA3 + STAT6 double, holding on the simulated response and the highlighted causal path
each time.

*Narration:* "MMC reconstructs a runnable model of the Th2 circuit from the newest CD4 T-cell
atlas, and it is not a black box: you can intervene on it. Knock down GATA3 and the model
simulates the cytokines falling and traces the causal path. Knock down STAT6 and the effect
travels through GATA3. Knock down both, and the model shows the second knockdown adds almost
nothing, because STAT6 acts through GATA3; the additive baseline sums the two and over-predicts,
while the mechanistic model composes them. A model you can interrogate is the capability the
black-box predictors do not provide."

---

**2. Grounded, not predictive (1:05-2:00)**

*On screen:* the "Held-out evaluation" view: the STK11 proposal and rationale, then the
scaled-corpus figure and metrics.

*Narration:* "The loop also proposes new hypotheses. On the cytokine module it proposes STK11 as
a chemokine regulator, and the STK11 edges are real: an edge-ablation control flags them
predictively necessary, like textbook edges. But across a powered corpus of 76 distinct novel
hypotheses over nine runs, none produces a model that beats a simple linear baseline held-out.
Grounded is not the same as predictive, and the property that separates them is held-out
advantage over a baseline, not the interpretability checks."

---

**3. Held to the baseline (2:00-2:35)**

*On screen:* the held-out result table declining the STK11 model (0.18 versus 0.45).

*Narration:* "So the engine holds its own proposal to the baseline and declines it. STK11 is a
real effect, but the model built from it does not predict better than the baseline, and the
engine reports it as proposed, not certified: a model that declines its own good-looking
hypotheses on the evidence."

---

**4. Close (2:35-3:00)**

*On screen:* the "Limit map" view.

*Narration:* "Single-perturbation prediction does not beat simple baselines, on this atlas and
in the field's own benchmarks, and the combinatorial regime does not rescue it. That is the case
for an instrument that can be interrogated and that audits itself against a baseline, rather than
one more prediction number."

---

Numbers: the intervention responses are from `demo/demo_interventions.json`; the corpus figures
from `scripts/engineer_behavior_scaled.py`; the held-out result from the cytokine evaluation
(0.18 versus 0.45). No prediction win and no disease discovery is claimed.
