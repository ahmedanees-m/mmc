"""The pre-registered metrics."""
import numpy as np

from mmc.eval import metrics


def test_sign_accuracy_is_taken_over_de_genes_only():
    pred = np.array([1.0, -1.0, 0.2])
    obs = np.array([2.0, -3.0, 0.0])       # only the first two are differentially expressed
    assert metrics.sign_accuracy(pred, obs, 1.0) == 1.0


def test_pearson_perfect_and_degenerate():
    x = np.array([1.0, 2.0, 3.0])
    assert abs(metrics.pearson(x, 2 * x) - 1.0) < 1e-9
    assert metrics.pearson(x, np.zeros(3)) is None      # a constant vector has no correlation


def test_de_overlap_top_k():
    pred = np.array([3.0, 0.1, -2.0, 0.0])
    obs = np.array([2.5, 0.0, -3.0, 0.2])
    ov = metrics.de_overlap(pred, obs, 2)               # top-2 by magnitude are {0, 2} in both
    assert ov["precision_at_k"] == 1.0 and ov["jaccard"] == 1.0


def test_bootstrap_ci_brackets_high_correlation():
    x = np.arange(20.0)
    ci = metrics.bootstrap_ci(metrics.pearson, x, x, n=200, seed=1)
    assert ci["lo"] > 0.9


def test_score_set_reports_pooled_and_per_perturbation():
    genes = ["G1", "G2", "G3"]
    pred = {"P": {"G1": 1.0, "G2": -1.0, "G3": 0.0}}
    obs = {"P": {"G1": 2.0, "G2": -2.0, "G3": 0.1}}
    r = metrics.score_set(pred, obs, genes, threshold=1.0)
    assert r["n_perturbations"] == 1
    assert r["pooled"]["sign_accuracy"] == 1.0
    assert "P" in r["per_perturbation"]
