# MMC

[![CI](https://github.com/ahmedanees-m/mmc/actions/workflows/ci.yml/badge.svg)](https://github.com/ahmedanees-m/mmc/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/ahmedanees-m/mmc/branch/main/graph/badge.svg)](https://codecov.io/gh/ahmedanees-m/mmc)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

MMC (Mechanistic Model Compiler) is a discovery loop for gene-regulatory circuits. It reads a
genome-scale perturbation atlas, proposes interpretable and runnable circuit models through a
language-model reasoning step, fits them to the measured perturbation responses, and evaluates
each model against simple baselines on held-out data. The reasoning step sets the structure and
logic of a model; a deterministic optimizer sets its numerical parameters.

The repository contains the loop, a deterministic evaluation harness, pre-registrations, the
result artifacts, and a Streamlit presentation of the results.

## Contents

- [Overview](#overview)
- [How it works](#how-it-works)
- [Results](#results)
- [Design decisions](#design-decisions)
- [Installation](#installation)
- [Usage](#usage)
- [Reproducibility](#reproducibility)
- [Repository layout](#repository-layout)
- [Data and references](#data-and-references)
- [License](#license)

## Overview

**What it is.** MMC pairs a language-model reasoning step with a deterministic modeling and
evaluation harness. The reasoning step proposes an executable model of a gene-regulatory module
as signed edges and bounded logic gates, reads the systematic residuals after fitting, and
revises the structure. The harness compiles the model to a steady state, fits its parameters,
and measures held-out predictive accuracy against simple baselines.

**What it produces.** For a module, MMC produces an interpretable circuit model, a
leave-one-perturbation-out evaluation of that model against mean and linear baselines, an
auditable trace of every proposal and revision, and a placement of the module on a limit map of
where mechanistic and AI models beat simple baselines. The reconstructed circuit is
interrogable: an intervention is simulated mechanistically, composing perturbations in a way an
additive baseline cannot (see the demo).

**Key finding.** The failure mode of AI-driven mechanistic discovery observed here is not
fabrication. The loop's novel hypotheses are individually supported and pass edge-level
necessity and ablation checks like established biology, but across a powered corpus of 76
distinct novel hypotheses none composes into a model that beats a linear baseline held-out
(Wilson 95% CI [0, 4.8%]). Interpretability checks, including necessity and ablation, are not
the safeguard; held-out advantage over a strong baseline is.

**Scope.** Results are on the Zhu 2025 genome-scale CD4+ T-cell Perturb-seq atlas, the modules
listed below, and single-knockdown steady-state readouts, with the combinatorial regime
evaluated on the Norman 2019 K562 dataset. No prediction win and no disease discovery is claimed.

## How it works

The loop separates two responsibilities: the reasoning step sets structure and logic; the
optimizer sets magnitudes.

```
Module (genes and differential expression, training state)
  -> Reasoning step: propose an executable model (signed edges, bounded logic gates);
     read the structural residuals; revise the structure, not the parameters
  -> Structural backend: compile to a steady-state fixed point; apply do(x=0) knockdown
     and do(x=high) activation interventions; fit parameters by gradient descent (JAX)
  -> Residual reader and ensemble: iterate to a plateau; freeze an ensemble within a loss margin
  -> Held-out evaluation: leave-one-perturbation-out against mean, linear, and zero baselines
     on ACC_DEG and DE-overlap, with bootstrap confidence intervals
```

**Grammar.** A model is a set of signed edges and, per target gene, a bounded sum-of-products of
sigmoid gates (an interpretable disjunctive normal form). A single additive term is the monotone
default; a product of gates encodes AND; a sum of terms encodes OR; a negative weight encodes
NOT. Bounds (at most three terms per target, at most three regulators per term) keep the model
interpretable and identifiable.

**Backend.** Each gene's steady-state level is a bounded function of its regulators, solved as a
damped fixed point rather than by integrating dynamics, because the data are interventional
steady states. A knockdown is a `do(x = 0)` clamp and an activation is a `do(x = high)` clamp.
Parameters are fit by automatic differentiation through the vectorized fixed point.

**Evaluation.** The held-out evaluation refits the structure on all but one perturbation and
predicts the held-out perturbation, comparing the model against mean and linear baselines on
DE-overlap and sign accuracy over differentially expressed genes, with bootstrap confidence
intervals. An edge-ablation variant measures whether removing an edge lowers held-out accuracy.

## Results

**Circuit reconstruction.** On the Th2 / GATA3 module the loop reconstructs a coherent circuit
(mutual GATA3 / TBX21 antagonism, GATA3 amplification of IL4/IL5/IL13) at in-sample Pearson 0.93.

**Held-out evaluation.** On the cytokine-production module the loop proposes STK11 / LKB1 as a
chemokine regulator. An edge-ablation control confirms the STK11 edges are individually
supported: removing them lowers the model's held-out accuracy, comparable to textbook edges such
as GATA3 -> IL5 (`paper/GATE_DISCRIMINATION.md`, `paper/gate_discrimination.json`). The model
built from these edges nonetheless does not beat a linear baseline on held-out data (DE-overlap
0.18 versus 0.45, separated confidence intervals). Across the immune modules, edge-level support
does not imply predictive advantage over a simple baseline; the module-level held-out evaluation
is the check that distinguishes them.

**Scaled characterization.** Across 25 proposals and 76 distinct novel hypotheses (110
instances) accumulated over 9 runs and two conditions, none is in a model that beats a linear
baseline held-out (validated rate 0 of 76, Wilson 95% CI [0, 4.8%]); 0 of 9 module-conditions
beat linear (CI [0, 30%]) (`scripts/engineer_behavior_scaled.py`,
`paper/engineer_behavior_scaled.png`). A reasoning-versus-search comparison shows the reasoning
step is not equivalent to random structure search: the proposed structure fits in-sample better
than a random structure of equal edge count on the informative modules (cytokine 0.19 versus
0.00, Th2 0.23 versus 0.08 Pearson; mean 0.20 versus 0.07), and not on the redundant TCR cascade
(0.39 versus 0.46), where structure carries little information.

The module-level held-out evaluation is a discriminator, not a uniform rejecter: on synthetic
ground truth it certifies a true structure (held-out DE-overlap 0.90 versus a mean baseline's
0.17) and ranks a fully-connected over-connected structure below it (0.63)
(`scripts/gate_synthetic_control.py`, `paper/gate_synthetic_control.json`).

**Limit map.** MMC places each module on a map of where mechanistic and AI models beat simple
baselines on held-out data (`paper/LIMIT_MAP.md`, `paper/mmc_limit_map.png`). On single-knockdown
steady-state data the mechanistic model shows no held-out advantage over simple baselines in any
measured regime, whether or not it fits in-sample. The map is consistent with independent
benchmarks: Ahlmann-Eltze et al. 2025 for single perturbations, and combinatorial benchmarks for
double perturbations.

**Combinatorial regime.** The remaining source of advantage is non-additivity, which
single-knockdown data cannot exercise. A compose test on the Norman 2019 K562 CRISPRa dataset
fits a structural model on single perturbations and predicts held-out double perturbations
(`scripts/norman_pseudobulk.py`, `scripts/norman_epistasis.py`, `paper/NORMAN_RESULT.md`). On the
non-additive pairs the model does not beat an additive baseline (held-out DE-overlap 0.35 versus
0.37), because the pair-specific interaction is not identifiable from single-perturbation
marginals.

**Principal finding.** The failure mode observed here is not fabrication. The loop's novel
hypotheses are individually supported and pass edge-level necessity and ablation checks like
established biology, and the reasoning step produces better-fitting structure than random search
on the informative modules. What the hypotheses do not do is compose into a model that beats a
simple baseline on held-out data: across 76 distinct novel hypotheses none does (Wilson 95% CI
[0, 4.8%]). Support at the level of individual edges, including necessity and ablation, does not
imply predictive advantage over a baseline; only held-out predictive advantage over a strong
baseline establishes that. This distinguishes a supported marginal effect from a model whose
predictions can be relied upon.

## Design decisions

- **Steady state over dynamics.** Fitting ordinary-differential-equation dynamics to
  steady-state knockdown data is under-identified. The structural steady-state backend fits the
  same data directly and is verified by a self-consistency check that recovers a known model at
  Pearson 1.0 (`scripts/structural_selfcheck.py`).
- **Held-out over in-sample.** In-sample fit does not imply held-out prediction; a strong
  in-sample fit on the Th2 module (0.93) does not hold under leave-one-perturbation-out. All
  reported model-versus-baseline comparisons are held-out.
- **Pre-registration.** Success rules are fixed before each run in `prereg/`.
- **Data-driven module selection.** Modules and the cytokine-module gene set are selected by a
  precondition and baseline-beatability screen (`scripts/baseline_screen.py`,
  `scripts/cytokine_module.py`), not by hand.

## Installation

```
pip install -e .            # core
pip install -e ".[dev]"     # tests and linting
pip install -e ".[diffrax]" # JAX backend for structural fitting
```

Requires Python 3.11 or later.

## Usage

The Streamlit presentation runs on precomputed results and needs no data or key:

```
pip install streamlit
streamlit run demo/app.py
```

Running the loop against the atlas requires the differential-expression store and a model API
key:

```
export MMC_ZHU_STORE=/path/to/zhu/store   # perturbation.parquet, gene.parquet, de.parquet
export ANTHROPIC_API_KEY=...              # read from the environment, never committed
python scripts/run_discovery.py <module> <condition> <iterations>
```

The reasoning model is configurable with `MMC_MODEL` and defaults to `claude-opus-4-8`.

## Reproducibility

- The test suite is offline and requires no store, key, or GPU: `pytest -q` (42 tests).
- Pre-registrations are committed before their runs: `prereg/PREREG_discovery.md` (cytokine
  discovery) and `prereg/PREREG_norman.md` (combinatorial). `prereg/PREREG.md` is the superseded
  transfer pre-registration, retained as a record.
- Result artifacts and figures are in `paper/`. The held-out evaluation is `mmc/eval/holdout.py`.
- The atlas is queried, not re-ingested, via `MMC_ZHU_STORE`. No keys or tokens are stored in the
  repository.

## Repository layout

```
mmc/            core library
  grammar/      model specification (signed edges, bounded logic gates)
  compile/      structural steady-state backend and interventions
  fit/          parameter fitting and the fit-versus-structure diagnostic
  loop/         propose, read residuals, repair, ensemble
  eval/         held-out evaluation, metrics
  baselines/    mean, linear, consensus baselines
  data/         module extraction, precondition, splits
  shared/       store accessor, model client, receipts, trace
scripts/        pipelines and analyses
tests/          offline unit tests
paper/          result artifacts and figures
prereg/         pre-registrations
demo/           Streamlit presentation
```

## Data and references

- Zhu et al. 2025, genome-scale Perturb-seq of primary human CD4+ T cells (bioRxiv
  10.64898/2025.12.23.696273).
- Norman et al. 2019, combinatorial CRISPRa Perturb-seq in K562 (GEO GSE133344).
- Moonen et al. 2026, variant, enhancer, and gene mapping in CD4+ T cells (bioRxiv
  10.64898/2026.03.09.710372).
- Ahlmann-Eltze, Huber, and Anders 2025, Nature Methods (10.1038/s41592-025-02772-6).

## License

MIT. See [LICENSE](LICENSE).
