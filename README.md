# mmc

[![CI](https://github.com/ahmedanees-m/mmc/actions/workflows/ci.yml/badge.svg)](https://github.com/ahmedanees-m/mmc/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24363795.svg)](https://doi.org/10.5281/zenodo.24363795)

Code for a structure-source comparison on genome-scale single-cell perturbation data. Each
source of regulatory structure, whether proposed by a language model, taken from canonical
immunology, inferred by a network-inference method, or found by direct search, is compiled into
one bounded sum-of-products grammar, fitted with one objective, and scored on one
leave-one-perturbation-out protocol against a regularized linear map.

On this atlas no structure source outperforms that linear map. A search permitted to select on
the folds it is scored on reaches a median 88% of its own score on data whose
perturbation-to-response pairing has been destroyed, which places most of its apparent advantage
with selection rather than with recovered structure.

## Installation

```bash
pip install -e .
```

Python 3.11 or later. `Dockerfile` builds the CPU image, `Dockerfile.gpu` adds the CUDA runtime
for the JAX backend, and `Dockerfile.gears` builds the environment for the graph-prior
comparator.

## Data

None of the source data is redistributed here.

- **CD4+ T cell CRISPRi Perturb-seq**, Zhu et al., *Cell* 189 (2026) 6988-7012,
  doi:10.1016/j.cell.2026.08.002. The input read here is the authors' deposited
  differential-expression summary statistics, `GWCD4i.DE_stats.h5ad`, from their CZ Virtual
  Cell Models deposit at
  https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq. The
  underlying sequencing data are in GEO under accession GSE314342, which does not carry the
  summary statistics.
- **K562 CRISPRa Perturb-seq**, Norman et al. 2019, GEO GSE133344. Used for the combinatorial
  regime, pseudobulked by `scripts/norman_pseudobulk.py`.
- **CollecTRI** regulons, OmniPath export, for the regulon-derived modules and the
  annotation-agreement analysis.

Analyses on the CD4+ atlas read per-module response matrices by one of two routes.

| Environment variable | Input |
|---|---|
| `MMC_MODULE_MATRICES` | a directory of exported matrices, one `<module>__<condition>.npz` per module and cell state, each holding the effect-size matrix, the FDR matrix and the gene and perturbation labels |
| `MMC_ZHU_STORE` | a directory holding `perturbation.parquet`, `gene.parquet` and `de.parquet`, to the schema documented at the top of `mmc/shared/store.py` |

`MMC_MODULE_MATRICES` is checked first and the store is the fallback, so a deposited copy of
the matrices is enough to re-run any CD4+ analysis without the atlas.
`scripts/export_module_matrices.py` writes that directory from the store.

Those matrices, the pseudobulked K562 matrix and every result file underlying the paper's
tables and figures are archived at Zenodo under
[10.5281/zenodo.24363795](https://doi.org/10.5281/zenodo.24363795).

`data/module_definitions.json` carries every candidate module, its gene set and its screen
outcome, so module assembly can be inspected without either input.

## Running an analysis

Each entry point writes JSON and takes `--help`.

| Entry point | What it runs | Needs |
|---|---|---|
| `export_module_matrices.py` | exports per-module matrices from the store | store |
| `assemble_modules.py` | assembles candidate modules and applies the power screen | matrices |
| `run_comparator.py` | scores every structure source on one module under the shared protocol | matrices |
| `run_permutation_control.py` | the same search on data with the perturbation pairing destroyed | result JSON |
| `run_nested_search.py` | per-fold selection against the linear map, under three fold pairings | result JSON |
| `pool_nested_folds.py` | pools the outer folds and pairs at fold level | result JSON |
| `summarize_search_budget.py` | counts modules the search clears, on the seed maximum and median | result JSON |
| `compute_diagnostics.py` | identifiability diagnostics per module and cell state | matrices |
| `calibrate_generators.py` | sweeps simulated generators against the observed diagnostics | matrices |
| `residual_decomposition.py` | rescores every source on the perturbation-specific residual | matrices |
| `run_offset_arm.py` | rescores a recorded structure with a fitted per-gene offset added | matrices, result JSON |
| `summarize_offset_arms.py` | pools the offset arms across modules | result JSON |
| `rank_truncation.py` | compares a rank-truncated ridge against the full ridge | matrices |
| `summarize_truncation.py` | summarises the truncation results against two thresholds | result JSON |
| `run_proposer_arms.py` | the proposer arms, differing in what the proposal step is shown | matrices, model key |
| `run_proposer_panel.py` | the same proposal task across a panel of models | matrices, model key |
| `run_gears.py` | the graph-prior comparator against additive baselines | GEARS data |
| `norman_pseudobulk.py` | pseudobulks the combinatorial atlas | GSE133344 |
| `norman_combinatorial.py` | the combinatorial regime under the field-standard split | pseudobulk |

The six marked **result JSON** alone read result files only. They refit nothing and need no
input data, so each can be re-run against a copy of the results.

## Repository layout

| Path | Contents |
|---|---|
| `mmc/grammar/` | `model_spec.py`, the bounded sum-of-products grammar and its serialised form |
| `mmc/compile/` | `structural.py` and `structural_jax.py` (fixed-point solve and its JAX form), `perturb.py` (the clamp operator), `to_ode.py`, `simulate.py` |
| `mmc/fit/` | `fit_structural.py` (multi-start L-BFGS-B over the compiled loss), `fit_params.py`, `diagnose.py` |
| `mmc/data/` | `module_data.py` (the two input routes), `module_extract.py`, `module_generator.py` (candidate assembly and the power screen), `textbook.py` (canonical structures), `precondition.py` |
| `mmc/eval/` | `holdout.py` (the fold protocol and DE-overlap), `compare.py` (the paired statistic), `oracle_search.py` (the structure search), `random_null.py`, `identifiability.py`, `annotation.py`, `model_classes.py` |
| `mmc/baselines/` | `linear.py`, the ridge map and its reduced-rank form |
| `mmc/compare/` | `scp_infer_adapter.py`, the network-inference sources |
| `mmc/loop/` | `propose.py`, `repair.py`, `providers.py`, `anonymize.py`, `ensemble.py`, `residuals.py`, `run.py`, the proposal path |
| `mmc/shared/` | `store.py` (store schema and queries), `llm.py`, `trace.py` |
| `scripts/` | the entry points above |
| `tests/` | the offline test suite; no store, model key or GPU required |
| `data/` | `module_definitions.json` |
| `.github/workflows/` | lint and test on push |

## Citation

The archived code, derived inputs and result files are deposited at Zenodo under the concept
DOI [10.5281/zenodo.24363795](https://doi.org/10.5281/zenodo.24363795), which always resolves
to the current version.

```bibtex
@dataset{mmc,
  author    = {Mahaboob Ali, Anees Ahmed and Delhibabu, Radhakrishnan and Nelson, Everette Jacob Remington},
  title     = {Code and data from: Structure search inflates the apparent performance of interpretable gene regulatory models on genetic perturbation data},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24363795},
  url       = {https://doi.org/10.5281/zenodo.24363795}
}
```

## Licence

MIT. See [LICENSE](LICENSE).
