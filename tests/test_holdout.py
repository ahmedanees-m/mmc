"""Held-out (leave-one-perturbation-out) evaluation engine.

The engine is backend-agnostic: it takes fit and predict callables. These tests use trivial
mean-predicting callables and a small synthetic module, exercising the metrics, the LOO loop,
the baselines, edge ablation, and seed stability.
"""
import numpy as np

from mmc.eval.holdout import (
    ModuleData,
    acc_deg,
    de_overlap,
    edge_ablation_holdout,
    loo_evaluate,
    print_report,
    seed_stability,
)
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term


def test_acc_deg_sign_accuracy():
    pred = np.array([1.0, -1.0, 1.0])
    obs = np.array([1.0, -1.0, -1.0])
    de = np.array([True, True, True])
    assert acc_deg(pred, obs, de) == 2 / 3


def test_acc_deg_no_de_is_nan():
    assert np.isnan(acc_deg(np.array([1.0]), np.array([0.0]), np.array([False])))


def test_de_overlap_perfect_and_empty():
    pred = np.array([0.0, 0.0, 5.0])
    obs = np.array([0.0, 0.0, 2.0])
    de = np.array([False, False, True])
    assert de_overlap(pred, obs, de) == 1.0
    assert np.isnan(de_overlap(pred, obs, np.array([False, False, False])))


def _module():
    genes = ["A", "B", "C"]
    perts = ["A", "B"]
    observed = np.array([[0.0, -1.0, 2.0], [1.0, 0.0, -3.0]])
    de_mask = np.array([[False, True, True], [True, False, True]])
    spec = ModelSpec(
        genes=genes,
        edges=[Edge(regulator="A", target="B", sign=1),
               Edge(regulator="A", target="C", sign=1)],
        rules={"B": Rule(terms=[Term(regulators=["A"])]),
               "C": Rule(terms=[Term(regulators=["A"])])},
    )
    return ModuleData(genes, perts, observed, de_mask, spec)


def _backend():
    def fit_fn(spec, train_perts, train_obs, seed=0):
        return {"mean": np.asarray(train_obs).mean(axis=0)}

    def predict_fn(model, pert):
        return model["mean"]
    return fit_fn, predict_fn


def test_loo_evaluate_structure():
    mod = _module()
    fit_fn, predict_fn = _backend()
    rep = loo_evaluate(mod, fit_fn, predict_fn, n_boot=50, seed=0)
    for method in ("model", "mean", "zero"):
        assert method in rep
        assert len(rep[method]["acc_deg"]) == 3        # (mean, lo, hi)
        assert len(rep[method]["de_overlap"]) == 3
    print_report(rep)                                   # smoke: renders without error


def test_loo_evaluate_with_linear_baseline():
    mod = _module()
    fit_fn, predict_fn = _backend()

    def linear_fn(train_perts, train_obs, held):
        return np.zeros(len(mod.genes))
    rep = loo_evaluate(mod, fit_fn, predict_fn, linear_fn=linear_fn, n_boot=50)
    assert "linear" in rep


def test_edge_ablation_returns_rows():
    mod = _module()
    fit_fn, predict_fn = _backend()
    rows = edge_ablation_holdout(mod, fit_fn, predict_fn)
    assert len(rows) == len(mod.spec.edges)
    for r in rows:
        assert {"edge", "sign", "drop"} <= set(r)


def test_seed_stability():
    mod = _module()
    fit_fn, predict_fn = _backend()
    out = seed_stability(mod, fit_fn, predict_fn, seeds=(0, 1), edge=("A", "B"))
    assert len(out["holdout_acc_by_seed"]) == 2
    assert "ablation_drop_by_seed" in out
