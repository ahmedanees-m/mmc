# A field guide to when mechanistic and AI modeling beats simple baselines

*Validated on the Zhu 2025 CD4+ T-cell Perturb-seq atlas and consistent with the field's
own benchmarks. Scope: CD4+ T cells, these modules, single-knockdown steady-state readouts,
with the combinatorial regime now measured on independent data (Norman K562).*

## The confusion this resolves

The field keeps asking whether mechanistic and foundation models beat simple baselines at
perturbation prediction, and getting contradictory answers. Ahlmann-Eltze, Huber, and Anders
(Nature Methods 2025, 10.1038/s41592-025-02772-6) show deep and foundation models do not beat
a linear or mean baseline on single-perturbation prediction. Combinatorial benchmarks
(PerturBench; the Norman 2019 genetic-interaction analysis) show a fitted-additive baseline
beats foundation models on most double perturbations, because most dual effects are close to
linear. Yet mechanistic models clearly help somewhere. The map below reconciles this: models
beat baselines only in a specific, characterizable regime, and nowhere else.

## The map (measured on the atlas)

Two axes decide whether a mechanistic or AI model beats a simple baseline on held-out data:

1. **In-sample fittability**: can the model class even represent the responses? Low
   fittability is a data ceiling (the specific signal is swamped by the bulk response); it is
   necessary but not sufficient.
2. **Held-out advantage over the best simple baseline (mean, linear)**: the only thing that
   matters. This is what the map plots (`mmc_limit_map.png`).

The measured regimes, held-out DE-overlap (leave-one-perturbation-out), model versus the best
simple baseline:

| Regime | Module (this atlas) | Fittable? | Beats baseline held-out? | Why |
|---|---|---|---|---|
| Redundant cascade | TCR signalosome | yes | no (mean unbeatable) | perturbations 0.91 correlated; the mean already is the answer, no headroom |
| Weak specific signal | CD4 TF network | no (~0.2) | no | the specific TF-to-TF signal is swamped by the bulk response; a data ceiling |
| Strong fit, no prediction | Th2/GATA3; cytokine module | yes (in-sample 0.93) | **no** (model 0.18 < linear 0.45, separated CIs) | fitting is not predicting; capacity does not compose to held-out generalization |
| Non-additive (combinatorial) | Norman K562 (independent, measured) | yes (represents epistasis) | **no** (model 0.35 vs additive 0.37 on held-out doubles, non-additive pairs) | the interaction is unidentifiable from singles: a logic-gate model fit on the singles predicts doubles no better than additive even where non-additivity is large (0.77) |

The figure encodes provenance honestly: filled markers are measured held-out DE-overlap on
this atlas, open markers are other-metric or qualitative placements, and the star marks the
Norman combinatorial regime. That regime is now measured directly (`NORMAN_RESULT.md`): a
structural model fit on the singles does not beat additive on the held-out doubles, so even
the one regime with headroom in principle yields no held-out advantage from singles.

## The synthesis

A mechanistic or AI model beats a simple baseline on held-out prediction only where all three
hold, and each is often violated:

- **(a) The specific signal is strong**: not swamped by the bulk, non-specific response. On
  single-knockdown steady-state data this frequently fails (the CD4 TF network), and even
  when it holds (Th2/GATA3, the cytokine module) the model still fits without predicting.
- **(b) There is non-redundancy**: headroom over the mean. Functionally redundant modules
  (the TCR signalosome) leave none: the mean already is the answer.
- **(c) There is non-additivity**: the distinctive feature of mechanism (logic gates, the
  sum-of-products) is only needed for epistasis. The field's combinatorial benchmarks show
  most dual effects are close to additive, so the headroom is thin; and on the Norman atlas we
  measured that even where non-additivity is large, a logic-gate model fit on the singles does
  not beat additive on the held-out doubles, because the interaction is unidentifiable from the
  single-perturbation marginals.

The negative that follows: on single-knockdown steady-state data, mechanism has no held-out
advantage over simple baselines in any regime, whether or not it can fit. The only genuine
headroom is non-additivity; single-knockdown data cannot access it, and on combinatorial data
a model informed only by the singles still cannot predict it (the Norman compose test,
`NORMAN_RESULT.md`), because the interaction is unidentifiable from the single-perturbation
marginals. Mechanism fits but does not predict, in both the single-perturbation and the
combinatorial regime. This is why the field's answers conflict: models are usually evaluated
where nothing can win, and where they can fit they still do not predict.

## Why this is useful regardless of any positive

The map tells a target team, before they trust a model on a new module, whether that module is
even in a regime where mechanism can help: is the signal specific, is there redundancy, is
there non-additivity. It converts a contradictory literature into a single reproducible
boundary, validated on the newest atlas. A tool you can trust because it declares where it
cannot be trusted.
