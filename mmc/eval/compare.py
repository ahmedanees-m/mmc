"""Structure-source comparison on one shared set of held-out folds (PREREG_v4 section 2).

Every structure source in the Step 1 comparator, whether its edges came from the
reasoning step, from a textbook, from a causal-discovery algorithm, from a random
draw, or from a search with access to the answer, is compiled into the same
structural backend, fit with the same optimizer budget, and scored on the same
leave-one-perturbation-out folds. The only thing that varies is where the edges
came from, which is what makes the comparison a statement about the data rather
than about any one proposer.

The advantage statistic here is paired. `holdout.loo_evaluate` bootstraps each
method's fold scores independently and leaves the pairing on the table, so two
methods that track each other fold by fold can still show overlapping marginal
intervals. Perturbations differ enormously in how predictable they are, and that
between-fold variance is shared by every method, so differencing within a fold
removes it. PREREG_v4 section 1.2 fixes the paired form as the test: an advantage
counts only when the lower bound of the paired interval clears zero.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .holdout import ModuleData, acc_deg, de_overlap

METRICS = ("de_overlap", "acc_deg")
N_BOOT = 10_000

# The optimizer budget every structural source is held to. Fixed here rather than
# passed per call so that no source can be given a larger search by accident.
FIT_STARTS = 4
FIT_MAX_ITER = 250


@dataclass
class SourceResult:
    """Per-fold scores for one structure source, plus what produced them."""

    name: str
    predictions: np.ndarray                      # (n_perts, n_genes)
    scores: dict[str, np.ndarray]                # metric -> (n_perts,)
    n_edges: int = 0
    meta: dict = field(default_factory=dict)


def fold_predictions(mod: ModuleData, source_fn) -> np.ndarray:
    """Leave-one-perturbation-out predictions from a source.

    source_fn(train_idx, held_idx) returns the predicted delta over mod.genes for
    the held-out perturbation, having seen only the training rows.
    """
    n = len(mod.perts)
    out = np.zeros_like(mod.observed)
    for i in range(n):
        train = [j for j in range(n) if j != i]
        out[i] = np.asarray(source_fn(train, i), dtype=float)
    return out


def score_predictions(mod: ModuleData, preds: np.ndarray) -> dict[str, np.ndarray]:
    """Per-fold metric values. NaN where the fold has no observed DE gene."""
    n = len(mod.perts)
    return {
        "de_overlap": np.array(
            [de_overlap(preds[i], mod.observed[i], mod.de_mask[i]) for i in range(n)], float
        ),
        "acc_deg": np.array(
            [acc_deg(preds[i], mod.observed[i], mod.de_mask[i]) for i in range(n)], float
        ),
    }


def bootstrap_mean(vals: np.ndarray, *, n_boot: int = N_BOOT, seed: int = 0) -> dict:
    """Marginal mean and percentile interval over folds, for the reported tables."""
    v = np.asarray(vals, float)
    v = v[~np.isnan(v)]
    if v.size == 0:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan"), "n_folds": 0}
    rng = np.random.default_rng(seed)
    boots = v[rng.integers(0, v.size, size=(n_boot, v.size))].mean(axis=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"mean": float(v.mean()), "lo": float(lo), "hi": float(hi), "n_folds": int(v.size)}


def paired_delta(source: np.ndarray, comparator: np.ndarray, *,
                 n_boot: int = N_BOOT, seed: int = 0) -> dict:
    """The PREREG_v4 section 1.2 advantage statistic.

    Mean over folds of the within-fold difference, with a paired bootstrap over
    folds. `advantage` is True only when the interval's lower bound exceeds zero.
    Folds where either score is undefined are dropped from both arms together, so
    the arms stay aligned.
    """
    a = np.asarray(source, float)
    b = np.asarray(comparator, float)
    if a.shape != b.shape:
        raise ValueError(f"fold arrays differ in shape: {a.shape} vs {b.shape}")
    ok = ~(np.isnan(a) | np.isnan(b))
    d = a[ok] - b[ok]
    if d.size == 0:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "n_folds": 0, "advantage": False}
    rng = np.random.default_rng(seed)
    boots = d[rng.integers(0, d.size, size=(n_boot, d.size))].mean(axis=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"delta": float(d.mean()), "lo": float(lo), "hi": float(hi),
            "n_folds": int(d.size), "advantage": bool(lo > 0)}


def benjamini_hochberg(pvals: list[float], q: float = 0.05) -> list[bool]:
    """BH step-up. Returns the reject flags in the input order (PREREG_v4 section 1.3)."""
    p = np.asarray(pvals, float)
    n = p.size
    if n == 0:
        return []
    order = np.argsort(p)
    ranked = p[order]
    thresh = q * (np.arange(1, n + 1) / n)
    passing = np.flatnonzero(ranked <= thresh)
    out = np.zeros(n, bool)
    if passing.size:
        out[order[: passing[-1] + 1]] = True
    return out.tolist()


def permutation_p(source: np.ndarray, comparator: np.ndarray, *,
                  n_perm: int = N_BOOT, seed: int = 0) -> float:
    """One-sided p for "source beats comparator", by sign-flipping the paired differences.

    Used only to feed the BH families; the reported claim rests on the paired
    interval, not on this p.
    """
    a = np.asarray(source, float)
    b = np.asarray(comparator, float)
    ok = ~(np.isnan(a) | np.isnan(b))
    d = a[ok] - b[ok]
    if d.size == 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    flips = rng.choice([-1.0, 1.0], size=(n_perm, d.size))
    null = (flips * d).mean(axis=1)
    return float(((null >= d.mean()).sum() + 1) / (n_perm + 1))


# ------------------------- source constructors -------------------------
def bind_structural(spec, mod: ModuleData, *, n_starts: int = FIT_STARTS,
                    max_iter: int = FIT_MAX_ITER, activation: bool = False):
    """S1, S2, S4, S5: refit the given structure per fold and predict the held-out one.

    `activation` switches the intervention operator from do(x = 0) to do(x = high)
    for the CRISPRa datasets; the knockdown clamp is the default for Zhu.
    """
    from ..compile import structural
    from ..fit import fit_structural

    genes = list(mod.genes)
    intervene = structural.activation if activation else structural.knockdown

    def fn(train_idx, held_idx):
        observed = {
            mod.perts[i]: {genes[j]: float(mod.observed[i, j])
                           for j in range(len(genes)) if genes[j] != mod.perts[i]}
            for i in train_idx
        }
        fits = fit_structural.multi_fit(spec, observed, n_starts=n_starts, max_iter=max_iter)
        return np.asarray(intervene(spec, fits[0]["params"], mod.perts[held_idx]))

    return fn


def bind_linear(mod: ModuleData, l2: float = 1.0):
    """S6: the regularized linear map, the bar every structure source has to clear."""
    from ..baselines import linear as linear_bl

    genes = list(mod.genes)

    def fn(train_idx, held_idx):
        td = {mod.perts[i]: {genes[j]: float(mod.observed[i, j]) for j in range(len(genes))}
              for i in train_idx}
        held = mod.perts[held_idx]
        pred = linear_bl.reconstruct(td, genes, [mod.perts[i] for i in train_idx],
                                     [held], l2=l2)[held]
        return np.array([pred.get(g, 0.0) for g in genes], float)

    return fn


def bind_mean(mod: ModuleData):
    """S7: the mean of the training perturbations (the Ahlmann-Eltze baseline)."""

    def fn(train_idx, held_idx):
        return mod.observed[list(train_idx)].mean(axis=0)

    return fn


def bind_zero(mod: ModuleData):
    """The floor: predict no change."""

    def fn(train_idx, held_idx):
        return np.zeros(len(mod.genes), float)

    return fn


def evaluate_source(mod: ModuleData, name: str, source_fn, *, n_edges: int = 0,
                    meta: dict | None = None) -> SourceResult:
    preds = fold_predictions(mod, source_fn)
    return SourceResult(name=name, predictions=preds,
                        scores=score_predictions(mod, preds),
                        n_edges=n_edges, meta=dict(meta or {}))


def comparator_table(results: dict[str, SourceResult], *, reference: str = "linear",
                     metric: str = "de_overlap", seed: int = 0) -> list[dict]:
    """One row per source: marginal mean with interval, and the paired delta vs reference."""
    ref = results[reference].scores[metric]
    rows = []
    for name, res in results.items():
        marg = bootstrap_mean(res.scores[metric], seed=seed)
        delta = paired_delta(res.scores[metric], ref, seed=seed)
        rows.append({
            "source": name,
            "n_edges": res.n_edges,
            f"{metric}_mean": round(marg["mean"], 4),
            f"{metric}_lo": round(marg["lo"], 4),
            f"{metric}_hi": round(marg["hi"], 4),
            "delta_vs_" + reference: round(delta["delta"], 4),
            "delta_lo": round(delta["lo"], 4),
            "delta_hi": round(delta["hi"], 4),
            "advantage": delta["advantage"],
            "n_folds": delta["n_folds"],
            **{k: v for k, v in res.meta.items()},
        })
    return sorted(rows, key=lambda r: -r[f"{metric}_mean"])


def jaccard_edges(spec_a, spec_b) -> float:
    """Unsigned edge-set Jaccard, for the structure-agreement matrix."""
    ea = {(e.regulator, e.target) for e in spec_a.edges}
    eb = {(e.regulator, e.target) for e in spec_b.edges}
    if not ea and not eb:
        return float("nan")
    return len(ea & eb) / len(ea | eb)
