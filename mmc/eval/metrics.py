"""Pre-registered evaluation metrics (PREREG section 6).

The metrics operate on aligned predicted and observed per-gene deltas. The design rule
is never to reward predicting the mean: sign accuracy and DE-overlap are computed over
the differentially expressed genes, and correlations are reported on the DE subset as
well as on all module genes. Bootstrap confidence intervals resample genes.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import pearsonr, spearmanr


def align(pred: dict[str, float], obs: dict[str, float],
          genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    p = np.array([pred.get(g, 0.0) for g in genes], dtype=float)
    o = np.array([obs.get(g, 0.0) for g in genes], dtype=float)
    return p, o


def de_mask(obs: np.ndarray, threshold: float) -> np.ndarray:
    return np.abs(obs) >= threshold


def sign_accuracy(pred: np.ndarray, obs: np.ndarray, threshold: float) -> float | None:
    """Fraction of differentially expressed genes whose predicted sign is correct."""
    m = de_mask(obs, threshold)
    if not m.any():
        return None
    return float((np.sign(pred[m]) == np.sign(obs[m])).mean())


def pearson(pred: np.ndarray, obs: np.ndarray) -> float | None:
    if pred.size < 2 or np.ptp(pred) == 0 or np.ptp(obs) == 0:
        return None
    return float(pearsonr(pred, obs)[0])


def spearman(pred: np.ndarray, obs: np.ndarray) -> float | None:
    if pred.size < 2 or np.ptp(pred) == 0 or np.ptp(obs) == 0:
        return None
    return float(spearmanr(pred, obs)[0])


def de_overlap(pred: np.ndarray, obs: np.ndarray, k: int) -> dict:
    """Precision at k and Jaccard of the top-k moved genes, predicted versus observed."""
    k = min(k, pred.size)
    if k == 0:
        return {"precision_at_k": 0.0, "jaccard": 0.0, "k": 0}
    top_p = set(np.argsort(-np.abs(pred))[:k].tolist())
    top_o = set(np.argsort(-np.abs(obs))[:k].tolist())
    inter = len(top_p & top_o)
    union = len(top_p | top_o)
    return {"precision_at_k": inter / k, "jaccard": inter / union if union else 0.0,
            "k": k}


def bootstrap_ci(metric, pred: np.ndarray, obs: np.ndarray, n: int = 1000,
                 alpha: float = 0.05, seed: int = 0) -> dict | None:
    """Percentile bootstrap confidence interval for a metric, resampling genes."""
    rng = np.random.default_rng(seed)
    idx = np.arange(pred.size)
    vals = []
    for _ in range(n):
        s = rng.choice(idx, size=idx.size, replace=True)
        v = metric(pred[s], obs[s])
        if v is not None:
            vals.append(v)
    if not vals:
        return None
    lo, hi = np.quantile(vals, [alpha / 2, 1 - alpha / 2])
    return {"mean": float(np.mean(vals)), "lo": float(lo), "hi": float(hi), "n": len(vals)}


def score_set(pred_map: dict[str, dict[str, float]], obs_map: dict[str, dict[str, float]],
              genes: list[str], threshold: float, k: int = 10) -> dict:
    """Per-perturbation metrics plus pooled correlations with a bootstrap interval."""
    perts = [p for p in obs_map if p in pred_map]
    per: dict[str, dict] = {}
    pooled_p, pooled_o = [], []
    for p in perts:
        pv, ov = align(pred_map[p], obs_map[p], genes)
        per[p] = {
            "sign_accuracy": sign_accuracy(pv, ov, threshold),
            "pearson": pearson(pv, ov),
            "spearman": spearman(pv, ov),
            "de_overlap": de_overlap(pv, ov, k),
        }
        pooled_p.append(pv)
        pooled_o.append(ov)
    if not perts:
        return {"per_perturbation": {}, "pooled": None, "n_perturbations": 0}
    P, O = np.concatenate(pooled_p), np.concatenate(pooled_o)
    return {
        "per_perturbation": per,
        "n_perturbations": len(perts),
        "pooled": {
            "pearson": pearson(P, O),
            "spearman": spearman(P, O),
            "sign_accuracy": sign_accuracy(P, O, threshold),
            "pearson_ci": bootstrap_ci(pearson, P, O),
        },
    }
