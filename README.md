# mmc

[![CI](https://github.com/ahmedanees-m/mmc/actions/workflows/ci.yml/badge.svg)](https://github.com/ahmedanees-m/mmc/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

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

- **CD4+ T cell CRISPRi Perturb-seq**, Zhu et al. 2025, doi:10.64898/2025.12.23.696273. The
  analysis queries a prebuilt DuckDB and Parquet store rather than re-ingesting the atlas. Set
  `MMC_ZHU_STORE` to a directory holding `perturbation.parquet`, `gene.parquet` and `de.parquet`;
  the schema each must follow is documented at the top of `mmc/shared/store.py`. Building that
  store from the deposited `GWCD4i.DE_stats.h5ad` is outside this repository.
- **K562 CRISPRa Perturb-seq**, Norman et al. 2019, GEO GSE133344. Used for the combinatorial
  regime, pseudobulked by `scripts/norman_pseudobulk.py`.
- **CollecTRI** regulons, OmniPath export, for the regulon-derived modules and the
  annotation-agreement analysis.

`data/module_definitions.json` carries every candidate module, its gene set and its screen
outcome, so module assembly can be inspected without the atlas.

## Running an analysis

Thirteen entry points, one per analysis. Each writes JSON and takes `--help`.

| Entry point | What it runs |
|---|---|
| `assemble_modules.py` | assembles candidate modules from the atlas and applies the power screen |
| `run_comparator.py` | scores every structure source on one module under the shared protocol |
| `run_permutation_control.py` | repeats the search on data with the perturbation pairing destroyed |
| `run_nested_search.py` | selects structure on inner folds and scores it on held-out outer folds |
| `compute_diagnostics.py` | identifiability diagnostics per module and cell state |
| `calibrate_generators.py` | sweeps simulated generators against the observed diagnostics |
| `residual_decomposition.py` | rescores every source on the perturbation-specific residual |
| `rank_truncation.py` | compares a rank-truncated ridge against the full ridge |
| `run_proposer_arms.py` | the proposer arms differing in what the proposal step is shown |
| `run_proposer_panel.py` | the same proposal task across a panel of models |
| `run_gears.py` | the graph-prior comparator against additive baselines |
| `norman_pseudobulk.py` | pseudobulks the combinatorial atlas |
| `norman_combinatorial.py` | the combinatorial regime under the field-standard split |

## Repository layout

| Path | Contents |
|---|---|
| `mmc/` | the library: grammar, compilation, fitting, the proposal loop, evaluation, baselines |
| `scripts/` | the thirteen entry points above |
| `tests/` | the offline test suite, no store, key or GPU required |
| `data/` | module definitions |
| `.github/workflows/` | lint and test on push |

## Citation

```bibtex
@software{mmc,
  author  = {Mahaboob Ali, Anees Ahmed},
  title   = {mmc: structure-source comparison on single-cell perturbation data},
  year    = {2026},
  url     = {https://github.com/ahmedanees-m/mmc},
  license = {MIT}
}
```

## Licence

MIT. See [LICENSE](LICENSE).
