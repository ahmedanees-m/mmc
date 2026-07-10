# MMC: an AI biological engineer that knows when it is wrong

MMC is a Claude-driven loop that reads a genome-scale immune-perturbation atlas,
autonomously proposes interpretable, runnable models of gene-regulatory circuits, tests
them against the real perturbation data, and **refuses to certify the hypotheses it cannot
stand behind**, catching its own plausible-but-wrong ideas.

The honest headline: on this data, mechanistic models **fit but do not predict** held-out
single-knockdown responses, and MMC says so. The contribution is not a prediction win. It
is a working, self-correcting mechanistic-hypothesis engine, plus a rigorous map of where
mechanism can and cannot be trusted.

> Scope: the Zhu 2025 genome-scale CD4+ T-cell Perturb-seq atlas, these modules, CD4+ T
> cells. No prediction win and no disease discovery is claimed.

## The demo (three beats)

```
pip install streamlit
streamlit run demo/app.py
```

1. **It builds.** Claude writes a runnable, interpretable model of the Th2/GATA3 axis
   (allergy, asthma, atopic disease) from the atlas, as a signed simulatable circuit, and
   fits it (in-sample Pearson 0.93).
2. **It catches itself.** Claude proposes STK11/LKB1 as a chemokine repressor, reasons
   about it from the knockdown residuals, and the held-out gate **refuses to certify it**
   because it does not predict (held-out DE-overlap 0.18 versus a linear baseline's 0.45,
   cleanly separated CIs).
3. **It knows its limits.** The limit-map (`paper/mmc_limit_map.png`): mechanism has no
   held-out advantage over simple baselines on single-knockdown data in any regime, and the
   tool declares its boundary.

## What the method contributes

- **Reasoning-guided structural repair.** Executable gene-network synthesis exists (SCNS,
  CellNOpt, RACIPE), but those search structure space syntactically. MMC reads a residual
  pattern (a knockdown that raises a target the model lowers) and infers a structural change
  (flip an edge sign, insert a repressor). That inference is the step a syntactic search
  cannot perform.
- **Self-correction, by construction, and measured.** A held-out prediction gate and an
  anti-theater discovery protocol refuse to certify a hypothesis that fits but does not
  predict. On the cytokine module the loop proposed STK11 as a chemokine repressor with a
  coherent, data-grounded rationale; the gate rejected it. This was not a one-off: across
  every hypothesis the loop proposed on these modules, nine coherently-argued proposals
  yielded twenty-one novel edges and zero that survived held-out validation, a 100 percent
  gate catch rate. Plausibility did not track prediction, and the gate supplied the
  calibration the rationale lacked (`paper/engineer_behavior.png`). An AI that knows when it
  is wrong.
- **The limit-map.** A rigorous, mechanistically-explained boundary of where mechanism beats
  correlation and where it does not, resolving a live field confusion (why do mechanistic
  and foundation models keep failing to beat simple baselines? they are usually tested where
  nothing can, and fitting is not predicting).
- **An auditable model-evolution trace.** Every model version, structural edit, and
  rationale, including falsified hypotheses, is logged. The reasoning is an artifact.

## The iteration arc (what was wrestled with)

The project is an honest sequence of reframes, each forced by the data:

1. **Prediction-transfer framing (v1/v2).** Judge MMC on predicting unseen perturbations and
   cross-state transfer, decomposed into a conserved core and rewiring. Built the full
   pipeline: leakage-safe two-tier splits, mean/linear/consensus baselines, a fit-vs-
   structure gate, pre-registration.
2. **First honest negative.** On the TCR signalosome the model lost to the mean and
   persistence on every metric. A baseline-beatability screen showed why: the signalosome is
   the mean-baseline worst case (perturbations 0.91 correlated), so no method can win there.
3. **The model-class correction.** Fitting ODE *dynamics* to *steady-state* knockdown data
   is under-identified and slow (in-sample capped at Pearson ~0.30, unmoved by optimizer
   budget). Replaced the ODE backend with a structural steady-state model solved as a fixed
   point and fit by autodiff (JAX). A spec now fits in seconds; a self-consistency check
   recovers a known model at Pearson 1.0.
4. **In-sample is not prediction.** A strong in-sample fit on Th2 (0.93) collapsed on a
   proper held-out test. Built the held-out gate (leave-one-perturbation-out, model vs
   linear on DE-overlap with CIs) and a pleiotropy/novelty-filtered discovery module, and
   pre-registered the success rule before running.
5. **The decisive, clean negative.** On a powered, favorable-regime, rigorously-excluded
   cytokine module, the model predicts *worse* than a linear baseline held-out (DE-overlap
   0.18 vs 0.45, separated CIs). Decisive, not confounded: mechanism fits but does not
   predict single-knockdown steady-state data. The engineer's proposal (STK11) was logged as
   proposed-but-unvalidated — the discipline working.

Each step walked back an earlier overclaim rather than defend it. That is the craft.

## Architecture

The trust boundary is an actor and tool split: the reasoning step sets structure and logic
form; the optimizer sets magnitudes.

```
Module (genes + Zhu DE, training state)
  -> REASONING STEP (Claude): propose an executable structural model (signed edges + bounded
       logic gates), read the structural residuals, rewrite structure, not parameters
  -> STRUCTURAL BACKEND + GATE (deterministic): compile to a steady-state fixed point,
       do(x=0)/do(x=high) interventions, gradient-fit in JAX; a residual is structural only
       if the best fit gets a DE gene's direction wrong
  -> RESIDUAL READER + ENSEMBLE (deterministic): iterate to a plateau, freeze
  -> HELD-OUT GATE: leave-one-perturbation-out vs mean/linear/zero on ACC_DEG + DE-overlap
       with bootstrap CIs; certify only what predicts
```

## Reproducibility

- **The demo** runs on a clean clone with only Streamlit (it presents captured results; no
  API key, store, or GPU). `streamlit run demo/app.py`.
- **Pre-registrations** are committed before their runs: `prereg/PREREG.md` (transfer),
  `prereg/PREREG_discovery.md` (cytokine discovery), `prereg/PREREG_norman.md` (epistasis).
- **The limit-map** figure is in `paper/`. **The held-out gate** is `mmc/eval/holdout.py`.
  The **module builder** and its Zhu-exclusion are `scripts/cytokine_module.py`.
- The atlas is queried, not re-ingested, via `MMC_ZHU_STORE`. No keys or tokens are in the
  repository. `pytest -q` runs the offline suite.

## Reuse from VERDICT

The modeling loop is new. The supporting layer is carried over and lives in `mmc/shared/`:
the Zhu DE store accessor (`store.py`), the model client with structured output and refusal
handling (`llm.py`), and the receipt and evolution-trace pattern (`receipts.py`, `trace.py`).
