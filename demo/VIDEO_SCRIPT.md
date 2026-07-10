# MMC: 3-minute demo video script

Record the voiceover and the screen capture separately; hold on the self-correction. Cut to
3:00 or under. Screen-capture each beat from `demo/app.py`. Never lead with the negative; the
hero is autonomy plus the centerpiece finding in Beat 2 (the AI's novel hypotheses are
grounded, pass the interpretability checks, and still do not compose into a model that beats a
baseline; only held-out predictive advantage catches that) plus honest limits.

---

**Beat 1 · It builds (0:00-1:00)**

*On screen:* the loop proposing, simulating, revising, then the rendered Th2/GATA3 circuit
(Beat 1 in the app).

*Voiceover:* "This is Claude reading the newest immune-disease atlas and writing a runnable,
interpretable model of the Th2 circuit that drives allergic and autoimmune disease: signed
edges, logic, a simulator, autonomously. Not a black box. Something a target team can read
and test."

---

**Beat 2 · It holds itself to the baseline (1:00-2:05).** Hold here.

*On screen:* Claude's STK11 proposal and its rationale, then the held-out gate refusing to
certify the model, then the engineer-behavior figure (grounded edges, model below baseline).

*Voiceover:* "Then it proposes a new hypothesis, STK11 as a chemokine repressor, with a
coherent mechanistic story it reads straight from the knockdown data. And here is the honest
part: the STK11 edges are real. An edge-ablation control flags them predictively necessary,
exactly like textbook edges. But when Claude assembles them into a mechanistic model and asks
whether that model beats a simple linear baseline on data it never saw, the answer is no,
0.18 against 0.45. So the engine reports a real effect and refuses to certify the model as a
discovery. Plausible, grounded mechanism is not predictive advantage, and it knows the
difference. An AI scientist that will not overclaim even its own grounded hypotheses, the
thing the field is most afraid it cannot do."

---

**Beat 3 · It knows its limits (2:05-2:45)**

*On screen:* the field-level limit-map (`paper/mmc_limit_map.png`, Beat 3 in the app).

*Voiceover:* "And it maps exactly where mechanism can be trusted and where it can't, a
boundary that matches the field's own benchmarks and resolves why these models keep failing
to beat simple baselines: they are usually tested where nothing can win, and fitting is not
predicting."

---

**Close (2:45-3:00)**

*On screen:* the target-team framing (the app footer).

*Voiceover:* "Interpretable, testable mechanistic hypotheses from the newest atlas, and an AI
that won't send a target team chasing the wrong one. Trust by construction."

---

The Beat 2 numbers are from `paper/gate_discrimination.json` (the STK11 edges are grounded)
and the cytokine held-out gate (0.18 versus 0.45). Framing rule: no prediction win, no disease
discovery; STK11 is shown grounded-but-not-a-certified-model, and you say exactly that.
