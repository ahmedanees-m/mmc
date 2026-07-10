# Demo walkthrough (3 minutes)

A narration and screen-capture guide for a short walkthrough of `demo/app.py`. Record the
narration and the screen capture separately and cut to 3:00 or under. Each section maps to one
view in the app.

---

**1. Circuit reconstruction (0:00-1:00)**

*On screen:* the loop proposing, simulating, and revising, then the reconstructed Th2 / GATA3
circuit (Circuit reconstruction view).

*Narration:* "MMC reads a genome-scale immune-perturbation atlas and produces a runnable,
interpretable model of the Th2 circuit, which is associated with allergic and autoimmune
disease: signed edges, logic, and a simulator. The structure is proposed by the reasoning step
and the magnitudes are set by the optimizer, so the model is inspectable and testable rather
than a black box."

---

**2. Held-out evaluation (1:00-2:05)**

*On screen:* the STK11 proposal and its rationale, the held-out result table, and the
engineer-behavior figure.

*Narration:* "On the cytokine-production module the loop proposes a new hypothesis, STK11 as a
chemokine repressor, derived from the knockdown residuals. The STK11 edges are individually
supported: an edge-ablation control flags them as predictively necessary, comparable to
textbook edges. The model built from them, however, does not beat a linear baseline on
held-out data, 0.18 against 0.45 with separated confidence intervals. The hypothesis is
therefore reported as proposed but not certified. Edge-level support is not the same as
predictive advantage over a baseline, and the held-out evaluation is what distinguishes them."

---

**3. Limit map (2:05-2:45)**

*On screen:* the limit map (`paper/mmc_limit_map.png`, Limit map view).

*Narration:* "MMC maps where a mechanistic or AI model does and does not beat simple baselines
on held-out data. The boundary is consistent with the field's benchmarks and explains why such
models frequently fail to beat baselines: they are commonly evaluated where no method has an
advantage, and in-sample fit does not imply held-out prediction."

---

**Close (2:45-3:00)**

*On screen:* the app footer.

*Narration:* "MMC produces interpretable, testable mechanistic hypotheses from a genome-scale
atlas and reports the regime within which its outputs are reliable."

---

Numbers in section 2 are from `paper/gate_discrimination.json` (edge-level support) and the
cytokine held-out evaluation (0.18 versus 0.45). No prediction win and no disease discovery is
claimed; STK11 is presented as proposed but not certified.
