# Demo

A Streamlit presentation of MMC results. It renders precomputed outputs, so it runs on a clean
clone with Streamlit alone and requires no API key, data store, or GPU.

```
pip install streamlit
streamlit run demo/app.py
```

Three views:

1. **Circuit reconstruction.** The reconstructed Th2 / GATA3 circuit (allergy, asthma, atopic
   disease) as a signed, simulatable graph, with the in-sample fit (Pearson 0.93).
2. **Held-out evaluation.** The STK11 / LKB1 chemokine-regulator hypothesis on the
   cytokine-production module, and the leave-one-perturbation-out result: the model predicts
   below a linear baseline (held-out DE-overlap 0.18 versus 0.45, separated CIs), so the
   hypothesis is reported as proposed but not certified.
3. **Limit map.** The regime boundary of where mechanistic and AI models beat simple baselines
   on held-out data, and where they do not.

Scope: the Zhu 2025 CD4+ T-cell Perturb-seq atlas and these modules. No prediction win and no
disease discovery is claimed.
