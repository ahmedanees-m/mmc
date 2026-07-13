# Demo

A Streamlit interface for MMC results, arranged as four tabs. It renders precomputed outputs, so
it runs on a clean clone with Streamlit alone and requires no API key, data store, or GPU.

```
pip install streamlit
streamlit run demo/app.py
```

The app opens at http://localhost:8501.

Four tabs:

1. **Discovery loop.** A deterministic replay of a captured run on the TCR signalosome core
   (Stim8hr). The reasoning step proposes the canonical activating cascade; the fit-versus-
   structure gate surfaces the structural residuals, which show that knockdowns raise their
   targets; the reasoning step revises the affected edges to repressive across three edits, and
   the training loss and structural-residual count fall (0.96 to 0.66, 10 to 6). Step through the
   iterations with "Next iteration". The trace is captured by `scripts/capture_loop_trace.py` and
   stored in `demo/loop_replay.json`; the reasoning text is shown verbatim.
2. **STK11 proposal.** The STK11 / LKB1 chemokine-regulator proposal and the held-out evaluation
   declining the model (DE-overlap 0.18 versus a linear baseline at 0.45, shown as a bar
   comparison). An edge-ablation control supports the individual edges, but the assembled model
   does not beat the baseline.
3. **Circuit.** The reconstructed Th2 / GATA3 circuit. Select an intervention (none, GATA3 KD,
   STAT6 KD, GATA3 + STAT6 KD, GATA3 activation) and the structural model simulates the per-gene
   response, highlights the causal path, and reports the predicted cytokines. Switching from GATA3
   KD to the GATA3 + STAT6 double holds IL5 at the mechanistic value while the additive baseline
   over-predicts it, because STAT6 acts through GATA3.
4. **Evidence.** The powered result (0 of 76 novel hypotheses beat the baseline held-out, Wilson
   95% CI [0, 4.8%]), a reasoning-versus-random comparison of in-sample fit, and the limit map.

The discovery-loop tab reflects one captured run. The loop is stochastic, so a fresh capture
reproduces the same qualitative arc (an activating cascade revised toward repression as the fit
improves) with different exact numbers.

Scope: the Zhu 2025 CD4+ T-cell Perturb-seq atlas and the modules shown.
