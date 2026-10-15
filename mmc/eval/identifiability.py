"""Identifiability diagnostics computed from the data alone (PREREG_v4 section 3).

These are the quantities that are meant to tell you, before any model is fit,
whether a module can support a structural model at all. They exist because "it did
not work" is not a reusable result, whereas "here is the measurable property of the
response matrix that predicts whether it can work" is.

All of them are computed on training rows only. The held-out rows are used in one
place, `leading_pc_variance_fraction`, and only as the thing being explained, never
as the thing fitted.

Two of the six are read off the S5 search trace rather than the raw data:
equivalence width and sign stability. Both are properties of the space of
near-optimal structures, which is exactly what non-identifiability is about, and
the search trace is where that space has been enumerated.
"""
from __future__ import annotations

import numpy as np


def _centered(matrix: np.ndarray) -> np.ndarray:
    m = np.asarray(matrix, float)
    return m - m.mean(axis=0, keepdims=True)


def effective_rank(responses: np.ndarray) -> float:
    """Participation ratio of the response matrix spectrum.

    (sum of eigenvalues)^2 / (sum of squared eigenvalues), on the covariance
    spectrum of the perturbation-by-gene response matrix. It runs from 1, meaning
    every perturbation moves the transcriptome along the same direction and nothing
    can tell them apart, up to the number of perturbations. This is the simplest
    candidate for the diagnostic that predicts the ceiling, and if it turns out to
    be the one that carries the signal then the result is easier to reuse.
    """
    x = _centered(responses)
    if x.shape[0] < 2:
        return float("nan")
    sv = np.linalg.svd(x, compute_uv=False)
    lam = sv ** 2
    total = lam.sum()
    if total <= 0:
        return float("nan")
    return float(total ** 2 / (lam ** 2).sum())


def leading_pc_variance_fraction(train: np.ndarray, held: np.ndarray) -> float:
    """Fraction of a held-out response's energy lying along the training PC1.

    This quantifies the stereotyped bulk response directly. When it is high, the
    mean-of-training baseline is already capturing most of what there is to
    capture, and a structural model has almost no room left to win in.
    """
    x = _centered(train)
    if x.shape[0] < 2:
        return float("nan")
    _, _, vt = np.linalg.svd(x, full_matrices=False)
    pc1 = vt[0]
    h = np.asarray(held, float) - np.asarray(train, float).mean(axis=0)
    denom = float(h @ h)
    if denom <= 0:
        return float("nan")
    return float((h @ pc1) ** 2 / denom)


def mean_leading_pc_fraction(responses: np.ndarray) -> float:
    """`leading_pc_variance_fraction` averaged over leave-one-out folds."""
    r = np.asarray(responses, float)
    n = r.shape[0]
    if n < 3:
        return float("nan")
    vals = []
    for i in range(n):
        train = r[[j for j in range(n) if j != i]]
        vals.append(leading_pc_variance_fraction(train, r[i]))
    vals = [v for v in vals if not np.isnan(v)]
    return float(np.mean(vals)) if vals else float("nan")


def perturbation_specific_ratio(responses: np.ndarray) -> float:
    """Perturbation-specific energy divided by shared energy.

    The response matrix is split into the common response, which is the same
    direction and magnitude for every perturbation, and what is left. A low ratio
    means the perturbations are mostly reproducing one shared program, which is the
    regime where no structure can distinguish them.
    """
    r = np.asarray(responses, float)
    if r.shape[0] < 2:
        return float("nan")
    shared = np.tile(r.mean(axis=0), (r.shape[0], 1))
    specific = r - shared
    shared_energy = float((shared ** 2).sum())
    if shared_energy <= 0:
        return float("nan")
    return float((specific ** 2).sum() / shared_energy)


def equivalence_width(train_loss, holdout, *, epsilon: float = 1.05) -> dict:
    """Spread of held-out score across structures that fit the training data equally well.

    Takes every structure within `epsilon` of the best training loss seen in the S5
    search and reports how widely their held-out scores range. A wide class at equal
    training fit is direct evidence that the data does not determine the structure,
    and it is the mechanism behind a ceiling: if many structures fit equally and
    predict differently, fitting well carries no information about predicting well.

    Epsilon is fixed at 1.05 in PREREG_v4 section 3 rather than chosen after seeing
    which value separates the modules.
    """
    loss = np.asarray(list(train_loss), float)
    hold = np.asarray(list(holdout), float)
    ok = ~(np.isnan(loss) | np.isnan(hold))
    loss, hold = loss[ok], hold[ok]
    if loss.size == 0:
        return {"n_in_class": 0, "width": float("nan")}
    best = float(loss.min())
    # Losses are non-negative; a multiplicative band is scale-free across modules.
    in_class = hold[loss <= best * epsilon] if best > 0 else hold[loss <= epsilon - 1.0]
    if in_class.size == 0:
        return {"n_in_class": 0, "width": float("nan")}
    return {
        "n_in_class": int(in_class.size),
        "n_evaluated": int(hold.size),
        "best_train_loss": best,
        "width": float(in_class.max() - in_class.min()),
        "sd": float(in_class.std(ddof=1)) if in_class.size > 1 else 0.0,
        "holdout_min": float(in_class.min()),
        "holdout_max": float(in_class.max()),
        "holdout_mean": float(in_class.mean()),
    }


def sign_stability(edge_sets, train_loss, *, epsilon: float = 1.05) -> dict:
    """Per-edge sign agreement across the near-optimal structures.

    The structural backend fixes an edge's sign in the structure and fits only a
    non-negative magnitude, so a sign cannot flip within one fit. The identifiable
    question is therefore whether the near-optimal structures agree on the sign of
    the edges they share. An edge that appears activating in half the equally
    well-fitting structures and repressing in the other half is an edge whose
    direction the data does not determine.
    """
    loss = np.asarray(list(train_loss), float)
    sets = list(edge_sets)
    if not sets or loss.size != len(sets):
        return {"n_edges": 0, "fraction_stable": float("nan")}
    best = float(np.nanmin(loss))
    band = [s for s, ln in zip(sets, loss)
            if not np.isnan(ln) and ln <= (best * epsilon if best > 0 else epsilon - 1.0)]
    if not band:
        return {"n_edges": 0, "fraction_stable": float("nan")}
    counts: dict[tuple[str, str], list[int]] = {}
    for s in band:
        for reg, tgt, sign in s:
            c = counts.setdefault((reg, tgt), [0, 0])
            c[0 if sign > 0 else 1] += 1
    if not counts:
        return {"n_edges": 0, "fraction_stable": float("nan")}
    stable, per_edge = 0, {}
    for (reg, tgt), (pos, neg) in counts.items():
        agree = max(pos, neg) / (pos + neg)
        per_edge[f"{reg}->{tgt}"] = round(agree, 3)
        if agree >= 0.9:
            stable += 1
    return {
        "n_edges": len(counts),
        "n_in_class": len(band),
        "fraction_stable": round(stable / len(counts), 3),
        "per_edge_agreement": per_edge,
    }


def effect_size_summary(responses: np.ndarray, de_mask: np.ndarray) -> dict:
    """Effect-size and knockdown-strength distribution, the input-side diagnostic.

    Weak interventions are the condition under which the causal-discovery
    literature reports failure, so this is recorded per module to place it against
    those benchmarks.
    """
    r = np.abs(np.asarray(responses, float))
    de = np.asarray(de_mask, bool)
    per_pert_de = de.sum(axis=1)
    de_eff = r[de]
    return {
        "n_perts": int(r.shape[0]),
        "n_genes": int(r.shape[1]),
        "n_de_entries": int(de.sum()),
        "de_per_pert_mean": float(per_pert_de.mean()),
        "de_per_pert_min": int(per_pert_de.min()),
        "de_per_pert_max": int(per_pert_de.max()),
        "abs_effect_median": float(np.median(r)) if r.size else float("nan"),
        "de_abs_effect_median": float(np.median(de_eff)) if de_eff.size else float("nan"),
        "de_abs_effect_p90": float(np.percentile(de_eff, 90)) if de_eff.size else float("nan"),
    }


def diagnostics(mod, *, trace: dict | None = None) -> dict:
    """All six diagnostics for one module. `trace` is the S5 search trace when available."""
    r = np.asarray(mod.observed, float)
    out = {
        "effective_rank": effective_rank(r),
        "effective_rank_normalised": effective_rank(r) / max(1, r.shape[0]),
        "leading_pc_fraction": mean_leading_pc_fraction(r),
        "perturbation_specific_ratio": perturbation_specific_ratio(r),
        "effect_sizes": effect_size_summary(r, mod.de_mask),
    }
    if trace:
        out["equivalence_width"] = equivalence_width(trace["train_loss"], trace["holdout"])
        if "edge_sets_raw" in trace:
            out["sign_stability"] = sign_stability(trace["edge_sets_raw"], trace["train_loss"])
    return out
