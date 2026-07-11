# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.8.1] - 2026-07-11

### Added
- Unit tests for the structural steady-state backend, the held-out evaluation engine, and the
  receipt, trace, and module-registration utilities (`tests/test_structural.py`,
  `tests/test_holdout.py`, `tests/test_shared.py`). Line coverage of the offline-testable core
  rises from 40 percent to 75 percent.

### Changed
- Coverage configuration targets the offline-testable core and excludes modules that require
  external resources unavailable in CI (the store, the model API, and the JAX backend).

## [0.8.0] - 2026-07-11

### Added
- Interactive circuit interrogation in the demo (`demo/app.py`, `demo/demo_interventions.json`,
  `scripts/demo_precompute.py`): the validated Th2 / GATA3 sub-circuit responds to interventions
  (knockdown, activation, and a double), simulating the per-gene response and tracing the causal
  path, and contrasting the mechanistic composition of a double knockdown with the additive
  baseline.
- Norman top-decile positive control (`scripts/norman_positive_control.py`,
  `paper/norman_positive_control.json`): the most strongly non-additive pairs show no held-out
  advantage over the additive baselines, and no per-pair case where the model beats both.

### Changed
- The demo leads with the interactive circuit; the video walkthrough is updated to match.
- The README overview elevates the key finding (grounded but not compositional) and the
  interrogability of the reconstructed circuit.

## [0.7.0] - 2026-07-10

### Added
- Scaled engineer-behavior characterization across modules, conditions, and repeated runs
  (`scripts/engineer_behavior_scaled.py`, `paper/engineer_behavior_scaled.json`,
  `paper/engineer_behavior_scaled.png`): the held-out validated rate and the module-condition
  beats-linear rate with Wilson confidence intervals, and a reasoning-versus-search comparison
  of the proposed structure against a random structure of equal edge count.

### Changed
- The engineer-behavior claim is reported with confidence intervals from the powered corpus
  (0 of 76 distinct novel hypotheses beat a linear baseline held-out, Wilson 95% CI [0, 4.8%])
  in the README, the abstract, and the demo.

## [0.6.0] - 2026-07-10

### Added
- Continuous integration (GitHub Actions) running Ruff and the test suite with coverage upload.
- Ruff configuration and coverage reporting via `pytest-cov`.
- `CHANGELOG.md`.

### Changed
- Documentation rewritten in a neutral technical register; `README.md` expanded into a full
  project description with installation, usage, results, and repository layout.
- Support layer (`mmc/shared/`) documented as native modules; external-repository references
  removed so the project is self-contained.
- `SUBMISSION.md` replaced by `ABSTRACT.md`.

## [0.5.0] - 2026-07-10

### Added
- Engineer-behavior characterization comparing proposal necessity against held-out predictive
  advantage (`scripts/engineer_behavior.py`, `paper/engineer_behavior.png`).
- Edge-ablation gate-discrimination control (`scripts/gate_discrimination.py`) and a synthetic
  discriminator control (`scripts/gate_synthetic_control.py`).
- Combinatorial compose test on the Norman 2019 K562 dataset: pseudobulk and epistasis
  evaluation (`scripts/norman_pseudobulk.py`, `scripts/norman_epistasis.py`,
  `paper/NORMAN_RESULT.md`).
- Limit-map field guide (`paper/LIMIT_MAP.md`) and result artifacts under `paper/`.
- Streamlit presentation of results (`demo/`).

### Changed
- Model-versus-baseline claims scoped to held-out predictive advantage; edge-level necessity
  distinguished from module-level advantage.

## [0.4.0] - 2026-07-10

### Added
- Cytokine-production module with a power precondition and program-based exclusion of known
  regulators (`scripts/cytokine_module.py`).
- CRISPRa activation operator `do(x = high)` with a toy verification
  (`scripts/activation_toycheck.py`).
- Pre-registered cytokine discovery run (`prereg/PREREG_discovery.md`).
- Limit-map figure (`paper/mmc_limit_map.png`).

### Fixed
- Proposal generation for targets with many candidate regulators (grammar bounds stated in the
  prompt).

## [0.3.0] - 2026-07-10

### Added
- Structural steady-state backend solved as a damped fixed point, fit by automatic
  differentiation (`mmc/compile/structural_jax.py`, `mmc/fit/fit_structural.py`).
- Leave-one-perturbation-out held-out evaluation with bootstrap confidence intervals and an
  edge-ablation variant (`mmc/eval/holdout.py`).

### Changed
- Replaced the ordinary-differential-equation dynamics backend with the structural steady-state
  backend on the fit path; parameter fitting no longer integrates dynamics.

## [0.2.0] - 2026-07-09

### Added
- Propose, read-residuals, repair, and ensemble loop (`mmc/loop/`).
- Two-tier split construction with a leakage audit, mean/linear/consensus baselines, DEG-level
  metrics, and negative and positive controls (`mmc/data/splits.py`, `mmc/baselines/`,
  `mmc/eval/`).
- Baseline-beatability screen for data-driven module selection (`scripts/baseline_screen.py`).

### Changed
- Fit objective weighted toward differentially expressed genes with a wrong-direction penalty.
- Fit-versus-structure diagnostic tightened to require agreement across restarts.

## [0.1.0] - 2026-07-09

### Added
- Model grammar: signed edges and bounded sum-of-products sigmoid gates
  (`mmc/grammar/model_spec.py`).
- Differential-expression store accessor and reasoning-step model client (`mmc/shared/`).
- Edge-level precondition test and data-driven module selection (`mmc/data/precondition.py`).
- Pre-registration template and an offline unit-test suite.
- Project packaging and Docker image.
