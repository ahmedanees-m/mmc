"""Step 5 and Step 6 evaluation: the in-context gate and the two-tier transfer test.

The frozen model predicts a perturbation by simulating its knockdown and reading the
steady-state delta. Step 5 scores the frozen model on the in-context held-out
perturbations against the mean and linear baselines, which is the go/no-go gate. Step 6
scores transfer to the test state: Tier A is strict transfer of the frozen model to
every test-state perturbation; Tier B scores a model that was adapted on the discovery
subset only, on the disjoint held-out subset. Nothing from the test state enters the
frozen model, and the Tier B held-out subset never enters the adaptation, which the
splits enforce by construction.

Scoring is kept free of any live loop: the Tier B adaptation is run by the caller and
its adapted model is passed in, so this module is deterministic and testable offline.
"""
from __future__ import annotations

import json

import numpy as np

from ..baselines import consensus as consensus_bl
from ..baselines import linear as linear_bl
from ..baselines import mean as mean_bl
from ..compile.perturb import knockdown
from ..compile.simulate import steady_state
from ..data import module_extract
from ..data.splits import build
from ..grammar.model_spec import ModelSpec
from . import metrics

DE_THRESHOLD = 0.5    # |log2 fold change| above which a gene is called differentially expressed


def load_frozen(path: str) -> tuple[ModelSpec, dict]:
    with open(path) as f:
        data = json.load(f)
    return spec_and_params(data["best"])


def spec_and_params(record: dict) -> tuple[ModelSpec, dict]:
    spec = ModelSpec.model_validate(record["spec"])
    p = record["params"]
    params = {"basal": np.array(p["basal"], float),
              "decay": np.array(p["decay"], float), "terms": p["terms"]}
    return spec, params


def model_predict(spec: ModelSpec, params: dict, perts: list[str],
                  genes: list[str]) -> dict[str, dict[str, float]]:
    """Predict each perturbation's per-gene delta by simulating its knockdown. The
    on-target self effect is dropped, matching the observed-delta convention."""
    gi = {g: i for i, g in enumerate(spec.genes)}
    wt = steady_state(spec, params)
    out: dict[str, dict[str, float]] = {}
    for p in perts:
        if p not in gi:
            continue
        d = knockdown(spec, params, p, wt=wt)
        out[p] = {g: float(d[gi[g]]) for g in genes if g in gi and g != p}
    return out


def _restrict(observed: dict, perts) -> dict:
    keep = set(perts)
    return {p: d for p, d in observed.items() if p in keep}


def _pooled_pearson(score: dict) -> float:
    pooled = score.get("pooled")
    if not pooled or pooled.get("pearson") is None:
        return float("nan")
    return float(pooled["pearson"])


def _ge(a: float, b: float) -> bool:
    """a is at least b, treating a missing baseline as not blocking."""
    return not np.isnan(a) and (np.isnan(b) or a >= b)


def insample(module: str, frozen_path: str, threshold: float = DE_THRESHOLD) -> dict:
    """Score the frozen model on the perturbations it was trained on. A low in-sample
    score means the fit did not converge, and any held-out result is uninterpretable
    until the optimizer budget is raised."""
    s = build(module)
    genes = module_extract.model_genes(module)
    spec, params = load_frozen(frozen_path)
    train_obs = _restrict(module_extract.observed_deltas(module, s.train_state), s.train_perts)
    pred = model_predict(spec, params, list(s.train_perts), genes)
    return {"stage": "insample", "scores": {"model": metrics.score_set(pred, train_obs, genes, threshold)}}


def incontext_gate(module: str, frozen_path: str, threshold: float = DE_THRESHOLD) -> dict:
    """Step 5: score the frozen model on the in-context held-out perturbations against
    the mean and linear baselines, in the training state."""
    s = build(module)
    genes = module_extract.model_genes(module)
    spec, params = load_frozen(frozen_path)

    train_obs = _restrict(module_extract.observed_deltas(module, s.train_state), s.train_perts)
    ic_obs = _restrict(module_extract.observed_deltas(module, s.train_state), s.incontext_heldout)

    preds = {
        "model": model_predict(spec, params, list(s.incontext_heldout), genes),
        "mean": mean_bl.predict(mean_bl.fit(train_obs, genes), list(s.incontext_heldout)),
        "linear": linear_bl.reconstruct(train_obs, genes, list(s.train_perts),
                                        list(s.incontext_heldout)),
    }
    scores = {k: metrics.score_set(v, ic_obs, genes, threshold) for k, v in preds.items()}
    pooled = {k: _pooled_pearson(v) for k, v in scores.items()}
    passed = _ge(pooled["model"], pooled["linear"]) and _ge(pooled["model"], pooled["mean"])
    return {"stage": "incontext_gate", "pooled_pearson": pooled,
            "passed_gate": bool(passed), "scores": scores}


def tier_a(module: str, frozen_path: str, threshold: float = DE_THRESHOLD) -> dict:
    """Step 6 Tier A: strict transfer of the frozen model to every test-state
    perturbation, against the mean and persistence baselines."""
    s = build(module)
    genes = module_extract.model_genes(module)
    spec, params = load_frozen(frozen_path)

    train_obs_all = module_extract.observed_deltas(module, s.train_state)
    train_obs = _restrict(train_obs_all, s.train_perts)
    test_obs = _restrict(module_extract.observed_deltas(module, s.test_state), s.tierA_perts)

    preds = {
        "model": model_predict(spec, params, list(s.tierA_perts), genes),
        "mean": mean_bl.predict(mean_bl.fit(train_obs, genes), list(s.tierA_perts)),
        # persistence uses the measured train-state effect of the same perturbation.
        "persistence": linear_bl.persistence(train_obs_all, list(s.tierA_perts)),
    }
    scores = {k: metrics.score_set(v, test_obs, genes, threshold) for k, v in preds.items()}
    pooled = {k: _pooled_pearson(v) for k, v in scores.items()}
    return {"stage": "tier_a", "pooled_pearson": pooled, "scores": scores,
            "beats_mean": _ge(pooled["model"], pooled["mean"]),
            "beats_persistence": _ge(pooled["model"], pooled["persistence"])}


def tier_b(module: str, adapted_spec: ModelSpec, adapted_params: dict,
           threshold: float = DE_THRESHOLD) -> dict:
    """Step 6 Tier B: score a test-state-adapted model on the held-out subset, against
    the consensus, mean, and persistence baselines. The adapted model was produced by the
    loop from the discovery subset only; only the discovery effects, which Tier B is
    allowed to see, enter the consensus baseline."""
    s = build(module)
    genes = module_extract.model_genes(module)

    train_obs_all = module_extract.observed_deltas(module, s.train_state)
    test_obs_all = module_extract.observed_deltas(module, s.test_state)
    held_obs = _restrict(test_obs_all, s.tierB_heldout)
    discovery_visible = _restrict(test_obs_all, s.tierB_discovery)

    preds = {
        "model": model_predict(adapted_spec, adapted_params, list(s.tierB_heldout), genes),
        "mean": mean_bl.predict(mean_bl.fit(_restrict(train_obs_all, s.train_perts), genes),
                                list(s.tierB_heldout)),
        "persistence": linear_bl.persistence(train_obs_all, list(s.tierB_heldout)),
        "consensus": consensus_bl.predict(
            train_obs_all, discovery_visible, genes, list(s.tierB_heldout),
            list(s.tierB_discovery), list(s.train_perts)),
    }
    scores = {k: metrics.score_set(v, held_obs, genes, threshold) for k, v in preds.items()}
    pooled = {k: _pooled_pearson(v) for k, v in scores.items()}
    return {"stage": "tier_b", "pooled_pearson": pooled, "scores": scores,
            "beats_consensus": _ge(pooled["model"], pooled["consensus"]),
            "beats_mean": _ge(pooled["model"], pooled["mean"])}
