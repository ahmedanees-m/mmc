# The combinatorial regime: a structural model fit on singles does not predict double perturbations

*Measured on the Norman 2019 CRISPRa genetic-interaction Perturb-seq atlas (K562, GSE133344).
Scope: this atlas, K562, gain-of-function; the capability claim only, no disease claim. This
note extends the limit-map ([LIMIT_MAP.md](LIMIT_MAP.md)) from the single-perturbation immune
regime to the combinatorial regime, and is pre-registered in
[prereg/PREREG_norman.md](../prereg/PREREG_norman.md).*

## The question

Track A found that on single-knockdown steady-state immune data a mechanistic model fits but
does not predict held-out perturbations better than a linear baseline. The one regime with
genuine headroom in principle is non-additivity (epistasis), which single-perturbation data
cannot exercise. Norman's combinatorial CRISPRa atlas can. The pre-registered test: does a
logic-gate structural model, informed only by the single perturbations, predict the held-out
double perturbations better than an additive baseline, specifically on the pairs where
non-additivity is real?

## Method

- **Pseudobulk.** The GSE133344 filtered matrix (33694 genes, 111668 cells, 362M nonzeros)
  was streamed and grouped by perturbation identity into 105 single-gene perturbations, 131
  double-gene perturbations, and a pooled non-targeting control, using the 91168 cells with
  good coverage and a single confident identity. Per-perturbation log2 fold-changes were taken
  against the pooled control. All 131 doubles have both of their singles present.
- **Pair selection (mechanical, pre-registered).** For each double (A, B), non-additivity was
  scored as the deviation of the observed double from the best fitted additive model
  c1*sA + c2*sB over the pair's DE genes (1 minus the additive fit's R-squared). The top
  tertile is the non-additive set (n=43, mean non-additivity 0.77), the bottom tertile the
  additive control (n=43, mean 0.40). The two tertiles separate cleanly. Operationalization
  note: the paper's precomputed genetic-interaction distance-correlation table is not shipped
  with the GEO matrix, so the same additive-fit deviation is recomputed from the pseudobulk;
  this is a stated deviation from the letter of the pre-registration.
- **Compose test (held-out).** The signed-logistic sum-of-products structural model (the mmc
  kernel, reduced to the bipartite A,B to readout topology) was informed only by the singles
  and the global response scale, then predicted the double via the activation operator with
  both inputs high, which saturates and is therefore non-additive. No double is seen at fit
  time. Baselines: fitted-additive (c1*sA + c2*sB fit to the observed double, the harder
  oracle bar) and mean-of-singles.
- **Gate.** DE-overlap (primary) and ACC_DEG on the held-out doubles, bootstrap confidence
  intervals across pairs, per set.

## Result

Held-out DE-overlap (precision-at-50), model versus fitted-additive and mean-of-singles:

| Set | n | non-additivity | model | fitted-additive | mean-of-singles |
|---|---|---|---|---|---|
| Non-additive | 43 | 0.77 | 0.349 [0.303, 0.396] | 0.376 [0.325, 0.430] | 0.370 [0.318, 0.422] |
| Additive control | 43 | 0.40 | 0.586 [0.541, 0.628] | 0.642 [0.605, 0.679] | 0.613 [0.571, 0.654] |

ACC_DEG follows the same pattern (non-additive: model 0.698, additive baselines 0.698 to 0.702;
additive control: model 0.862, additive baselines 0.862 to 0.870). The structural model does
not beat the additive baselines on either set and sits marginally below. This is the
pre-registered null.

A positive control on the most favorable subset confirms this. On the top decile of
non-additivity (13 pairs, deviation 0.79 to 0.99, the strongest genetic interactions),
held-out DE-overlap is model 0.19 versus fitted-additive 0.21 and mean-of-singles 0.20, and no
individual pair has the model beating both baselines (0 of 13). The most non-additive setting
is where a structural model would beat additive if it could anywhere, and it does not
(`scripts/norman_positive_control.py`, `paper/norman_positive_control.json`).

## What it means

The non-additive pairs are genuinely non-additive: even the fitted-additive oracle, which is
allowed to see the double, reaches only 0.38 DE-overlap on them (against 0.64 on the additive
controls). The headroom is real and large. But a model informed only by the singles cannot
reach it, because the pair-specific interaction is unidentifiable from single-perturbation
marginals. A saturating logic gate imposes one particular non-additive extrapolation; it is
not the true interaction, so it does not help and slightly hurts. The identifiability argument
is general for models that infer the double from the single-perturbation marginals alone, which
is what this compose test does. Methods that inject external pairwise priors (for example GEARS,
with its gene-gene graph) do modestly better on the non-additive subset, but even they do not
reliably beat a fitted-additive baseline (Ahlmann-Eltze, Huber, and Anders, Nature Methods
2025, report additive beating GEARS overall), so the headroom stays thin.

This closes the last "in principle" gap in the limit-map. The single genuine source of
headroom over simple baselines, non-additivity, is inaccessible not only to single-knockdown
data but also, from singles, in combinatorial data. Mechanism fits but does not predict, in
both the single-perturbation and the combinatorial regime. The negative is unifying, not a
one-off.

## Reproduce

`scripts/norman_pseudobulk.py` builds the pseudobulk from the GEO matrix (regenerated, not
committed; the raw atlas is 1.13 GB). `scripts/norman_epistasis.py` runs the pair selection,
compose test, and gate, and writes `paper/norman_result.json`. Both run in the mmc container.
