"""Algorithmic structure sources for the Step 1 S4 arm (PREREG_v4 section 2.2, amendment A3).

S4 was specified as seven interventional causal-discovery methods run through
`scp_infer`. Five of them cannot be run on this atlas, and the reason is a property of
the data rather than of the software.

`scp_infer` takes an AnnData of cells by genes with a per-cell perturbation label, and
GIES, DCDI-G, DCDI-DSF, Bicycle, AVICI and SDCD all estimate their scores from many
observations within each intervention environment. The Zhu store is not that. Its
manifest records it as derived from `GWCD4i.DE_stats.h5ad`, holding log fold change,
adjusted p value and z score per perturbation and gene: one summary row per
intervention, not a population of cells. GIES is the closest to runnable and still
fails, because its Gaussian score needs a non-singular covariance and one row per
environment gives n equal to p.

Two of the seven survive the reduction, and they are implemented here against the
response matrix directly:

  Mean Difference   the effect estimate itself, thresholded. This is what the store is.
  GRNBoost2         per-target gradient-boosted regression over the other genes, with
                    perturbations as samples. Faithful to the arboreto specification
                    that `scp_infer` wraps, implemented directly because at a 28 by 28
                    response matrix the wrapper's AnnData and Dask machinery adds a
                    large dependency tree and nothing else. Named accordingly wherever
                    it is reported.

The remaining five are not silently dropped. They are run on Norman in Step 4, where
per-cell counts exist, and the fact that a widely used comparator suite cannot be
applied to a published DE-summary atlas is reported as a finding in its own right.
"""
from __future__ import annotations

import numpy as np

from ..eval.holdout import ModuleData
from ..eval.oracle_search import EdgeKey

METHODS = ("mean_difference", "grnboost2")


def mean_difference_edges(mod: ModuleData, *, max_edges: int = 30,
                          require_de: bool = True) -> list[EdgeKey]:
    """Edges ranked by the size of the measured knockdown effect.

    The regulator is the perturbed gene, the target is the gene that moved, and the
    sign is the direction the target moved when the regulator was removed: a target
    that falls on knockdown is activated by it, so a negative log fold change gives a
    positive edge.
    """
    gi = {g: i for i, g in enumerate(mod.genes)}
    scored: list[tuple[float, EdgeKey]] = []
    for pi, pert in enumerate(mod.perts):
        for gene in mod.genes:
            if gene == pert:
                continue
            j = gi[gene]
            if require_de and not mod.de_mask[pi, j]:
                continue
            eff = float(mod.observed[pi, j])
            if eff == 0.0:
                continue
            scored.append((abs(eff), (pert, gene, 1 if eff < 0 else -1)))
    scored.sort(key=lambda kv: -kv[0])
    return _cap_in_degree([e for _, e in scored], max_edges)


def grnboost2_edges(mod: ModuleData, *, max_edges: int = 30, seed: int = 0,
                    n_estimators: int = 500, learning_rate: float = 0.01,
                    max_depth: int = 3, subsample: float = 0.9) -> list[EdgeKey]:
    """GRNBoost2 over the perturbation-response matrix.

    For each target gene, a gradient-boosting regressor is fit with the other genes'
    responses as features; the feature importances are the edge weights. The
    hyperparameters follow the arboreto GRNBoost2 defaults. Regulators are restricted
    to perturbed genes, since an unperturbed regulator cannot be exercised by the
    held-out folds. Edge signs are not produced by the importances, so they are taken
    from the mean-difference direction as the pre-registration specifies.
    """
    from sklearn.ensemble import GradientBoostingRegressor

    x = np.asarray(mod.observed, float)
    gi = {g: i for i, g in enumerate(mod.genes)}
    reg_idx = [gi[p] for p in mod.perts if p in gi]
    if len(reg_idx) < 2:
        return []

    weights: list[tuple[float, str, str]] = []
    for target in mod.genes:
        tj = gi[target]
        feats = [i for i in reg_idx if i != tj]
        if not feats:
            continue
        y = x[:, tj]
        if np.allclose(y, y[0]):
            continue
        model = GradientBoostingRegressor(
            n_estimators=n_estimators, learning_rate=learning_rate, max_depth=max_depth,
            subsample=subsample, random_state=seed)
        model.fit(x[:, feats], y)
        for f, imp in zip(feats, model.feature_importances_):
            if imp > 0:
                weights.append((float(imp), mod.genes[f], target))

    weights.sort(key=lambda kv: -kv[0])
    signs = _sign_lookup(mod)
    edges = [(r, t, signs.get((r, t), 1)) for _, r, t in weights]
    return _cap_in_degree(edges, max_edges)


def _sign_lookup(mod: ModuleData) -> dict[tuple[str, str], int]:
    """Edge signs from the mean-difference direction, per PREREG_v4 section 2.2."""
    gi = {g: i for i, g in enumerate(mod.genes)}
    pi = {p: i for i, p in enumerate(mod.perts)}
    out: dict[tuple[str, str], int] = {}
    for pert, i in pi.items():
        for gene, j in gi.items():
            if gene == pert:
                continue
            eff = float(mod.observed[i, j])
            if eff != 0.0:
                out[(pert, gene)] = 1 if eff < 0 else -1
    return out


def _cap_in_degree(edges: list[EdgeKey], max_edges: int,
                   max_in_degree: int = 3) -> list[EdgeKey]:
    """Threshold a ranked edge list into a grammar-legal binary edge set."""
    kept: list[EdgeKey] = []
    in_deg: dict[str, int] = {}
    seen: set[tuple[str, str]] = set()
    for reg, tgt, sign in edges:
        if len(kept) >= max_edges:
            break
        if (reg, tgt) in seen or in_deg.get(tgt, 0) >= max_in_degree:
            continue
        seen.add((reg, tgt))
        in_deg[tgt] = in_deg.get(tgt, 0) + 1
        kept.append((reg, tgt, sign))
    return kept


def algorithmic_edges(mod: ModuleData, method: str, **kw) -> list[EdgeKey]:
    if method == "mean_difference":
        return mean_difference_edges(mod, **kw)
    if method == "grnboost2":
        return grnboost2_edges(mod, **kw)
    raise ValueError(
        f"{method} needs per-cell interventional data, which this store does not hold; "
        f"runnable methods here are {METHODS}"
    )
