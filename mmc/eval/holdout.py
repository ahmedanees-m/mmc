"""Held-out (leave-one-perturbation-out) evaluation for MMC structural models.

Answers the question that gates Part 4: does the module's fit carry HELD-OUT
predictive signal, or is the high in-sample number just capacity? Compares the
structural model against the mean-of-training and zero baselines (and an optional
injected linear baseline) on the field's DEG-level metrics (ACC_DEG, DE-overlap),
with bootstrap CIs appropriate to small n.

Backend-agnostic by design: pass your own fit / predict callables. Bind points:
    fit_fn(spec, train_perts, train_obs, seed=0) -> fitted_model
    predict_fn(fitted_model, pert) -> np.ndarray  (predicted delta over `genes`)
Everything else (LOO loop, baselines, metrics, CIs, edge ablation) is here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ----------------------------- data container -----------------------------
@dataclass
class ModuleData:
    genes: list[str]                 # measured genes (columns of `observed`)
    perts: list[str]                 # perturbations = knocked-down genes (rows)
    observed: np.ndarray             # (n_perts, n_genes) observed delta (log_fc)
    de_mask: np.ndarray              # (n_perts, n_genes) bool: DE (significant) entries
    spec: object                     # the ModelSpec (structure) for the model path

    def __post_init__(self):
        self.observed = np.asarray(self.observed, float)
        self.de_mask = np.asarray(self.de_mask, bool)
        assert self.observed.shape == self.de_mask.shape
        assert self.observed.shape[0] == len(self.perts)
        assert self.observed.shape[1] == len(self.genes)


# ------------------------------- metrics ---------------------------------
def acc_deg(pred: np.ndarray, obs: np.ndarray, de: np.ndarray) -> float:
    """Sign accuracy over the observed DE genes (the field's ACC_DEG)."""
    idx = de & (np.abs(obs) > 0)
    if idx.sum() == 0:
        return np.nan
    return float((np.sign(pred[idx]) == np.sign(obs[idx])).mean())


def de_overlap(pred: np.ndarray, obs: np.ndarray, de: np.ndarray, k: int | None = None) -> float:
    """Jaccard of predicted top-|delta| genes vs observed DE set."""
    obs_de = set(np.flatnonzero(de).tolist())
    if not obs_de:
        return np.nan
    k = k or len(obs_de)
    pred_top = set(np.argsort(-np.abs(pred))[:k].tolist())
    return len(obs_de & pred_top) / len(obs_de | pred_top)


# ---------------------------- the LOO engine ------------------------------
def _predictions(mod: ModuleData, fit_fn, predict_fn, seed: int = 0,
                 spec=None) -> dict[str, np.ndarray]:
    """Return per-perturbation predictions for each method: (n_perts, n_genes)."""
    spec = spec if spec is not None else mod.spec
    n = len(mod.perts)
    out = {m: np.zeros_like(mod.observed) for m in ("model", "mean", "zero")}
    for i in range(n):
        train = [j for j in range(n) if j != i]
        train_obs = mod.observed[train]
        # model: refit structure on the other perturbations, predict the held-out one
        model = fit_fn(spec, [mod.perts[j] for j in train], train_obs, seed=seed)
        out["model"][i] = predict_fn(model, mod.perts[i])
        # baselines
        out["mean"][i] = train_obs.mean(axis=0)
        out["zero"][i] = 0.0
    return out


def _score(mod: ModuleData, preds: dict[str, np.ndarray]) -> dict[str, dict[str, np.ndarray]]:
    """Per-perturbation metric arrays for each method."""
    res = {m: {"acc_deg": [], "de_overlap": []} for m in preds}
    for m, P in preds.items():
        for i in range(len(mod.perts)):
            res[m]["acc_deg"].append(acc_deg(P[i], mod.observed[i], mod.de_mask[i]))
            res[m]["de_overlap"].append(de_overlap(P[i], mod.observed[i], mod.de_mask[i]))
        res[m] = {k: np.array(v, float) for k, v in res[m].items()}
    return res


def _boot_ci(vals: np.ndarray, n_boot: int, rng: np.random.Generator):
    v = vals[~np.isnan(vals)]
    if len(v) == 0:
        return (np.nan, np.nan, np.nan)
    means = [rng.choice(v, size=len(v), replace=True).mean() for _ in range(n_boot)]
    return float(np.mean(v)), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def loo_evaluate(mod: ModuleData, fit_fn, predict_fn, *, linear_fn=None,
                 n_boot: int = 2000, seed: int = 0) -> dict:
    """Leave-one-perturbation-out held-out evaluation. The Part-4 gate."""
    rng = np.random.default_rng(seed)
    preds = _predictions(mod, fit_fn, predict_fn, seed=seed)
    if linear_fn is not None:
        P = np.zeros_like(mod.observed)
        for i in range(len(mod.perts)):
            train = [j for j in range(len(mod.perts)) if j != i]
            P[i] = linear_fn([mod.perts[j] for j in train], mod.observed[train], mod.perts[i])
        preds["linear"] = P
    scored = _score(mod, preds)
    report = {}
    for m, mets in scored.items():
        report[m] = {met: _boot_ci(arr, n_boot, rng) for met, arr in mets.items()}
    report["_per_pert_acc_deg"] = {m: scored[m]["acc_deg"] for m in scored}
    return report


# ------------------- edge ablation on HELD-OUT (Part 4) -------------------
def _spec_without_edge(spec, regulator: str, target: str):
    """Return a copy of spec with one edge removed. Adapt to your ModelSpec API."""
    import copy
    s = copy.deepcopy(spec)
    s.edges = [e for e in s.edges if not (e.regulator == regulator and e.target == target)]
    # drop rule terms that referenced the removed regulator for this target
    if target in s.rules:
        for term in list(s.rules[target].terms):
            if regulator in term.regulators:
                term.regulators = [r for r in term.regulators if r != regulator]
        s.rules[target].terms = [t for t in s.rules[target].terms if t.regulators]
        if not s.rules[target].terms:
            del s.rules[target]
    return s


def edge_ablation_holdout(mod: ModuleData, fit_fn, predict_fn, *, seed: int = 0) -> list[dict]:
    """For each edge: held-out ACC_DEG drop when the edge is removed.

    A candidate edge is REQUIRED only if removing it lowers HELD-OUT ACC_DEG --
    not in-sample fit. This is the discipline Part 4 needs so a 'discovery' is a
    predictive necessity, not extra capacity.
    """
    base = loo_evaluate(mod, fit_fn, predict_fn, seed=seed)["model"]["acc_deg"][0]
    rows = []
    for e in mod.spec.edges:
        s = _spec_without_edge(mod.spec, e.regulator, e.target)
        preds = _predictions(mod, fit_fn, predict_fn, seed=seed, spec=s)
        abl = _score(mod, {"model": preds["model"]})["model"]["acc_deg"]
        abl_mean = float(np.nanmean(abl))
        rows.append({"edge": f"{e.regulator}->{e.target}", "sign": e.sign,
                     "holdout_acc_full": round(base, 3),
                     "holdout_acc_ablated": round(abl_mean, 3),
                     "drop": round(base - abl_mean, 3)})
    return sorted(rows, key=lambda r: -r["drop"])


def seed_stability(mod: ModuleData, fit_fn, predict_fn, seeds=(0, 1, 2, 3), *,
                   edge=None) -> dict:
    """Identifiability proxy: does the held-out ACC_DEG (and an edge's ablation
    drop) persist across refit seeds? A discovery must be seed-stable."""
    accs = [loo_evaluate(mod, fit_fn, predict_fn, seed=s)["model"]["acc_deg"][0] for s in seeds]
    out = {"holdout_acc_by_seed": [round(a, 3) for a in accs],
           "holdout_acc_std": round(float(np.std(accs)), 3)}
    if edge is not None:
        drops = []
        for s in seeds:
            s_minus = _spec_without_edge(mod.spec, *edge)
            full = loo_evaluate(mod, fit_fn, predict_fn, seed=s)["model"]["acc_deg"][0]
            abl = np.nanmean(_score(mod, {"model": _predictions(mod, fit_fn, predict_fn, seed=s, spec=s_minus)["model"]})["model"]["acc_deg"])
            drops.append(full - abl)
        out["edge"] = f"{edge[0]}->{edge[1]}"
        out["ablation_drop_by_seed"] = [round(d, 3) for d in drops]
        out["ablation_drop_std"] = round(float(np.std(drops)), 3)
    return out


def print_report(report: dict) -> None:
    print(f"{'method':<10}{'ACC_DEG (mean [95% CI])':<34}{'DE-overlap (mean [95% CI])'}")
    for m in ("model", "mean", "linear", "zero"):
        if m not in report:
            continue
        a = report[m]["acc_deg"]; d = report[m]["de_overlap"]
        print(f"{m:<10}{a[0]:.3f} [{a[1]:.3f}, {a[2]:.3f}]{'':<8}{d[0]:.3f} [{d[1]:.3f}, {d[2]:.3f}]")
