# MMC: Mechanistic Model Compiler

A mechanistic discovery loop that writes runnable, falsifiable models of T-cell
gene regulation and maps what is conserved versus rewired when the cell changes
state.

In an outer reasoning loop, a model writes a compact executable model of a
gene-regulatory subnetwork (a logic and ODE circuit as code). A deterministic
harness compiles and simulates it, reads the numerical failures against held-out
perturbation data, and returns the structural failures to the reasoning step,
which infers a mechanism and rewrites the model's structure. The loop iterates
toward a model that predicts unseen perturbations and, as the decisive test,
transfers across cell states: a model built on resting CD4+ T cells is used to
predict stimulated cells, and the circuit is decomposed into a conserved core and
a context-specific rewiring.

The governing rule: MMC is judged on prediction of unseen perturbations and
cross-state transfer, not on plausibility. A plausible mechanistic story is not
evidence. Every headline claim is a number on data the loop and the fitter never
saw.

## What the method contributes

- Reasoning-guided structural repair. Executable gene-network synthesis already
  exists (SCNS, CellNOpt, RACIPE), but those search structure space syntactically
  with constraint solvers or genetic algorithms. MMC's contribution is to read a
  residual pattern (for example, a knockdown that raises a target the model
  lowers) and infer a structural change (flip an edge sign, or insert an
  intermediate repressor). That inference is the step a syntactic search cannot
  perform.
- An auditable model-evolution trace. Every model version, every structural edit,
  and the rationale for it, including hypotheses that were later falsified, are
  logged. The reasoning is an artifact, not only the final model.
- Calibrated boundaries. The loop reports the conserved and rewired edges and
  returns "undetermined" where the data cannot distinguish structures.
- A reusable pattern. Executable-model synthesis, reasoning-guided structural
  repair, calibrated self-assessment, and cross-condition transfer validation form
  a loop that applies to any domain with a simulator and intervention data. The
  T-cell circuit is the first instance.

## Status

The decisive result for the project is Step 6 (cross-state transfer against
baselines, including Arc State), and that number can only come from the build. The
precondition that Step 6 depends on, whether there is a conserved core to transfer
and a rewiring to discover, has been measured on the Zhu atlas and is favorable.
Reproduced from the public Zhu summary table:

- The regulator-activity backbone is partially conserved across states:
  effect-breadth Spearman is 0.660 (Rest to Stim8hr), 0.602 (Rest to Stim48hr),
  and 0.656 (Stim8hr to Stim48hr) over 11,086 genes measured in all three states.
- Primary module, Th2 and GATA3. GATA3 is a conserved, active hub with strong
  training-state signal (it moves 887 downstream genes at Rest and has a clean
  on-target knockdown across all states); STAT6 provides a context-specific
  rewiring, moving 2, then 232, then 1127 downstream genes across Rest, Stim8hr,
  and Stim48hr.
- Second module, the TCR signalosome, run Stim8hr to Stim48hr. Its genes are
  effectively off at Rest (breadth 2 to 6), so a Rest-trained model would have no
  signal to learn from; the data confirms the direction.

This establishes that building is justified. It does not claim that MMC beats
Arc State; that remains the open question the build answers.

Reproduce with `python scripts/precondition_summary.py`. See `prereg/`.

## Architecture

The trust boundary is an actor and tool split: the reasoning step sets structure
and logic form, and the optimizer sets magnitudes. The same gate pattern that
separated computed values from asserted ones in VERDICT separates structure from
magnitude here.

```
Module spec (genes + Zhu DE summary, training state)
  -> REASONING STEP (synthesis + structural repair, logged)
       proposes an executable model (signed edges + logic/ODE rules),
       reads structural residuals, rewrites structure, not parameters
  -> COMPILER + SIMULATOR + FIT-VS-STRUCTURE GATE (deterministic)
       compile, simulate WT and knockdowns to steady state, multi-start fit;
       a residual is forwarded only if no seed can remove it
  -> RESIDUAL READER + ENSEMBLE TRACKER (deterministic)
       iterate to a plateau, then freeze the ensemble
  -> EVALUATE (two tiers): A strict transfer, B few-shot rewiring discovery
       against mean, linear, and Arc State baselines
       -> conserved/rewired map + bench-testable predictions
```

## Reuse from VERDICT

The modeling loop is new. The supporting layer is carried over from the existing
VERDICT build and lives in `mmc/shared/`:

- The Zhu DE store. MMC queries the existing DuckDB and Parquet store; it does not
  re-ingest the atlas. `mmc/shared/store.py` points at it via `MMC_ZHU_STORE`.
- The model client (structured output, adaptive thinking, refusal handling) in
  `mmc/shared/llm.py`.
- The receipt and trace pattern (provenance and the evolution log) in
  `mmc/shared/receipts.py` and `mmc/shared/trace.py`.
- The Moonen CRE-gene loader (orthogonal edge validation, Step 7) in
  `mmc/shared/moonen.py`.

What is new: `grammar/`, `compile/`, `fit/` (including the diagnostic gate),
`loop/` (propose, read, repair, ensemble), the baselines (including Arc State),
and the two-tier evaluation with the conserved and rewired decomposition.

## Win conditions (stated in advance)

- Strongest: matches or beats Arc State on cross-state transfer with an
  interpretable model and a validated conserved and rewired map.
- Strong: beats the linear baseline, is competitive with Arc State, and is the
  only method producing the interpretable decomposition, falsifiable predictions,
  and orthogonally corroborated edges.
- Honest negative: fails transfer, and the decomposition shows the module
  genuinely rewires. This is a reportable finding about context-specific
  regulation, demonstrated rather than only observed.

## Prior art

Executable gene-network synthesis exists (constraint solvers, 2015 onward); LLM
gene-network inference exists (static graphs, LLM4GRN); LLM executable-model
synthesis exists in ecology (LEMMA, the same method template); foundation models
now beat simple baselines on context transfer (Arc State, 2025); and LLM
mechanistic reasoning for perturbation exists (MechPert, SUMMER, SynthPert, AROMA),
producing neighborhoods or reasoning traces rather than runnable, structurally
repaired models. MMC is the first to use a reasoning model as the
biological-reasoning-guided structural-repair engine for executable, simulatable,
human-readable gene-regulatory models, validated by cross-state transfer and
decomposed into conserved versus rewired structure, on the newest genome-scale
primary human CD4+ T-cell atlas.

The premise is established: deep-learning perturbation predictors did not beat
simple linear or mean baselines (Ahlmann-Eltze, Huber, and Anders, Nature Methods
2025), and LLM mechanistic reasoning underperformed a gene-frequency baseline
(Plausibility Is Not Prediction, 2026). MMC answers both with a runnable,
falsifiable artifact.

## Layout

```
mmc/
  shared/     store  receipts  trace  llm  moonen     carried over from VERDICT
  data/       slice_zhu  module_extract  precondition  splits
  grammar/    model_spec                               schema + logic gates (implemented)
  compile/    to_ode  simulate  perturb                simulator core (implemented)
  fit/        fit_params  diagnose                     multi-start fit + fit-vs-structure gate
  loop/       propose  residuals  repair  ensemble
  baselines/  mean  linear  foundation (Arc State)  consensus
  eval/       metrics  controls  evaluate  conserved_rewired
prereg/       PREREG.md
tests/        grammar and simulator tests
```

## Quickstart

```bash
pip install -e ".[dev]"
pytest -q                                # grammar and simulator core tests
export MMC_ZHU_STORE=/path/to/zhu_store  # the existing Zhu Parquet store
export ANTHROPIC_API_KEY=...             # or an env file, chmod 600, never committed
```

## Datasets

Zhu et al. 2025, genome-scale Perturb-seq in primary human CD4+ T cells (bioRxiv
10.64898/2025.12.23.696273; CZI Virtual Cells Platform). Moonen et al. 2026,
variant to CRE to gene links across immune diseases (bioRxiv
10.64898/2026.03.09.710372), used for orthogonal edge validation. Arc State
(github.com/ArcInstitute/state), the foundation baseline.

## License

MIT. See `LICENSE`.
