# MMC demo

The three-beat demo of the honest AI biological engineer: Claude builds a runnable disease
circuit, catches its own plausible-but-wrong hypothesis, and declares its own limits.

The demo presents captured results, so it runs on a clean clone with only Streamlit; it
needs no API key, no data store, and no GPU.

```
pip install streamlit
streamlit run demo/app.py
```

The three beats:

1. **It builds**: Claude autonomously writes a runnable, interpretable model of the
   Th2 / GATA3 axis (allergy, asthma, atopic disease) from the atlas, rendered as a signed
   circuit, and fits it (in-sample Pearson 0.93).
2. **It catches itself**: Claude proposes STK11 / LKB1 as a metabolic repressor of
   chemokines, reasons about it from the knockdown residuals, and the held-out gate refuses
   to certify it because it does not predict (held-out DE-overlap 0.18 versus a linear
   baseline's 0.45, separated CIs).
3. **It knows its limits**: the limit-map: mechanism has no held-out advantage over simple
   baselines on single-knockdown data in any regime, and the tool says so.

No prediction win and no disease discovery is claimed. Every claim is scoped to the Zhu
2025 CD4+ T-cell atlas and these modules.
