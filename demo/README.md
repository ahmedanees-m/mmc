# Demo

A Streamlit presentation of MMC results, arranged as three screens for a short walkthrough. It
renders precomputed outputs, so it runs on a clean clone with Streamlit alone and requires no
API key, data store, or GPU.

```
pip install streamlit
streamlit run demo/app.py
```

Three tabs:

1. **Interrogate the circuit.** The reconstructed Th2 / GATA3 circuit. Choose an intervention
   (GATA3 KD, STAT6 KD, GATA3 + STAT6 KD, GATA3 activation) and the structural model simulates
   the response and traces the causal path. Switching from GATA3 KD to the GATA3 + STAT6 double
   holds IL5 at the mechanistic value while the additive baseline over-predicts it, because
   STAT6 acts through GATA3.
2. **The STK11 catch.** Claude's STK11 / LKB1 chemokine-regulator proposal, and the held-out
   gate declining the model (DE-overlap 0.18 versus a linear baseline's 0.45): grounded, not
   hallucinated, but not predictive.
3. **The evidence.** The powered result (0 of 76 novel hypotheses beat baseline held-out, Wilson
   95% CI [0, 4.8%]) and the limit map.

Scope: the Zhu 2025 CD4+ T-cell Perturb-seq atlas and these modules. No prediction win and no
disease discovery is claimed.
