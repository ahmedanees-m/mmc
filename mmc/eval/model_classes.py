"""Alternative model classes.

These separate two explanations for the result: that interpretable structure is too weak to
predict, or that the data does not support prediction by anything. The structural classes
all share the fixed-point solve and the `do(x = 0)` clamp, so each predicts identically zero
for a perturbation whose gene has no outgoing edges, and each therefore forfeits the shared
response that supplies most of the achievable signal.

The two classes here break that forfeit rather than varying the logic.

`offset` hands the structural model a fitted per-gene vector added to every perturbation's
prediction. That is exactly the component the model class cannot otherwise express, and
nothing more, so the comparison isolates it. If this closes the gap to the linear baseline,
the forfeit accounts for the gap.

`dense_linear` lifts the three-regulator cap on a signed linear structural model, testing
whether the cap alone binds, separately from the offset.

Both are thin wrappers over the existing fit path; a class that also changed the
optimizer would confound the comparison it is built to make.
"""
from __future__ import annotations

import numpy as np

from .compare import FIT_MAX_ITER, FIT_STARTS
from .holdout import ModuleData


def bind_structural_with_offset(spec, mod: ModuleData, *, n_starts: int = FIT_STARTS,
                                max_iter: int = FIT_MAX_ITER):
    """Structural prediction plus a per-gene offset fitted on the training folds.

    The offset is the mean training response, which is the mean baseline's entire
    prediction. Fitting it on the training rows of each fold keeps it out of fold: the
    held-out perturbation contributes nothing to the offset it is scored against.
    """
    from ..compile import structural
    from ..fit import fit_structural

    genes = list(mod.genes)

    def fn(train_idx, held_idx):
        observed = {
            mod.perts[i]: {genes[j]: float(mod.observed[i, j])
                           for j in range(len(genes)) if genes[j] != mod.perts[i]}
            for i in train_idx
        }
        fits = fit_structural.multi_fit(spec, observed, n_starts=n_starts,
                                        max_iter=max_iter)
        params = fits[0]["params"]
        structural_pred = np.asarray(structural.knockdown(spec, params,
                                                          mod.perts[held_idx]))
        # the component the structure cannot emit, estimated from training rows only
        train_struct = np.stack([
            np.asarray(structural.knockdown(spec, params, mod.perts[i]))
            for i in train_idx])
        offset = (mod.observed[list(train_idx)] - train_struct).mean(axis=0)
        return structural_pred + offset

    return fn


def bind_dense_linear(mod: ModuleData, *, l2: float = 1.0, max_regulators: int = 0):
    """A one-hot ridge over which gene was perturbed. **Not usable under leave-one-out.**

    Under leave-one-perturbation-out the held-out perturbation's indicator column never
    appears in the training rows, so its fitted coefficient is zero and the prediction is
    exactly zero for every gene. Measured, it reproduces the zero baseline exactly: 0.1220
    on CD4_lineage_TFs and 0.1252 on the cytokine module.

    Retained as the diagnostic that establishes this. An unbounded in-degree comparison
    has to be built as a structural model with the cap lifted, not as a regression on
    perturbation identity.

    `max_regulators` above zero keeps only that many largest-magnitude coefficients per
    target.
    """

    def fn(train_idx, held_idx):
        train = list(train_idx)
        n_genes = len(mod.genes)
        # one-hot design over which gene was perturbed
        x = np.zeros((len(train), n_genes))
        gi = {g: i for i, g in enumerate(mod.genes)}
        for r, i in enumerate(train):
            p = mod.perts[i]
            if p in gi:
                x[r, gi[p]] = 1.0
        y = mod.observed[train]
        a = x.T @ x + l2 * np.eye(n_genes)
        w = np.linalg.solve(a, x.T @ y)               # (n_genes perturbed, n_genes read)
        if max_regulators > 0:
            for col in range(w.shape[1]):
                order = np.argsort(-np.abs(w[:, col]))
                w[order[max_regulators:], col] = 0.0
        held = mod.perts[held_idx]
        vec = np.zeros(n_genes)
        if held in gi:
            vec[gi[held]] = 1.0
        return vec @ w

    return fn


def bind_mean_offset_only(mod: ModuleData):
    """The offset alone, with no structure. Identical to the mean baseline.

    Included so the table carries the decomposition explicitly: if
    structural-plus-offset matches this, the structural part contributes nothing beyond
    the offset, which is a sharper statement than comparing against linear alone.
    """

    def fn(train_idx, held_idx):
        return mod.observed[list(train_idx)].mean(axis=0)

    return fn
