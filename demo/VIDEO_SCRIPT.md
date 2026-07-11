# Demo walkthrough (3 minutes)

A narration and screen-capture guide for `demo/app.py`. Record the narration and the screen
capture separately and cut to 3:00 or under. Each section maps to one tab in the app.

---

**Tab 1: Interrogate the circuit (0:00-1:10)**

*On screen:* tab 1. Select GATA3 KD, then GATA3 + STAT6 KD. Hold on the two large IL5 cards: the
"MMC (mechanistic)" card stays at -5.0 while the "additive baseline" card moves to -5.4.

*Narration:* "MMC reconstructs a runnable model of the Th2 circuit from the newest CD4 T-cell
atlas, and it is not a black box: you can intervene on it. Knock down GATA3 and the model
simulates the cytokines falling and traces the causal path. Now add a second knockdown, STAT6.
The mechanistic IL5 response does not move, because STAT6 acts through GATA3, and once GATA3 is
off, knocking down STAT6 adds nothing. The additive baseline sums the two knockdowns and
over-predicts, reaching minus 5.4. Composing interventions is what a mechanistic model can do and
a baseline cannot. A model you can interrogate is the capability the black-box predictors do not
provide."

---

**Tab 2: The STK11 catch (1:10-2:10)**

*On screen:* tab 2. Claude's STK11 proposal on the left; the held-out gate on the right, the two
cards 0.18 versus 0.45 and the red NOT CERTIFIED.

*Narration:* "The loop also proposes new hypotheses. On the cytokine module it proposes STK11 as
a chemokine regulator, and the STK11 edges are real: an edge-ablation control flags them
predictively necessary, like textbook edges. But when Claude builds a model from them and asks
whether it beats a simple linear baseline on held-out data, the answer is no: 0.18 against 0.45.
The engine holds its own good-looking hypothesis to the baseline and declines it. Grounded is not
the same as predictive, and only the held-out comparison separates them."

---

**Tab 3: The evidence (2:10-2:45)**

*On screen:* tab 3. The two large cards, 0 of 76 and 0 of 9, and the limit map.

*Narration:* "This is not one example. Across a powered corpus, none of 76 distinct novel
hypotheses produces a model that beats a linear baseline held-out, and no module-condition does
either. The limit map places this in the field: single-perturbation prediction shows no advantage
over simple baselines in any measured regime, and the combinatorial regime does not rescue it."

---

**Close (2:45-3:00)**

*Narration:* "Prediction against simple baselines is stuck, in our results and in the field's own
benchmarks. That is the case for an instrument that can be interrogated and that audits itself
against a baseline, rather than one more prediction number."

---

Numbers: the intervention responses are from `demo/demo_interventions.json`; the corpus figures
from `scripts/engineer_behavior_scaled.py`; the held-out result from the cytokine evaluation (0.18
versus 0.45). No prediction win and no disease discovery is claimed.
