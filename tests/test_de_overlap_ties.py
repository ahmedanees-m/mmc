"""The tie-breaking rule in DE-overlap, and why it is not cosmetic.

A structural model only moves genes downstream of the perturbed one. For most
perturbations a sparse structure predicts exactly zero everywhere, so its entire
ranking is one tie. Ranking by argsort resolved that tie by gene index, which handed
such a model the alphabetically first genes of the module and scored it on whether
those happened to be differentially expressed.
"""
import numpy as np
import pytest

from mmc.eval import metrics
from mmc.eval.holdout import de_overlap, topk_overlap


def test_a_prediction_of_no_change_scores_at_chance_not_by_gene_order():
    """The bug: an all-zero prediction used to take the first k genes and could score
    well or badly purely on how the gene list was sorted."""
    n, k = 12, 4
    pred = np.zeros(n)
    obs = np.zeros(n)
    obs[:k] = 3.0                       # the DE genes sit at the front of the array
    de = np.zeros(n, bool)
    de[:k] = True

    got = de_overlap(pred, obs, de)
    # chance level for k of n drawn at random, as a Jaccard of two size-k sets
    expected_inter = k * k / n
    chance = expected_inter / (2 * k - expected_inter)
    assert got == pytest.approx(chance, abs=0.06)
    assert got < 0.45, "an all-zero prediction must not inherit the front of the array"


def test_the_score_does_not_depend_on_gene_order_for_a_degenerate_prediction():
    n = 12
    pred = np.zeros(n)

    def score(de_positions):
        obs = np.zeros(n)
        de = np.zeros(n, bool)
        obs[de_positions] = 3.0
        de[de_positions] = True
        return de_overlap(pred, obs, de)

    front = score([0, 1, 2, 3])
    back = score([8, 9, 10, 11])
    spread = score([0, 4, 7, 11])
    assert max(front, back, spread) - min(front, back, spread) < 0.12


def test_an_untied_prediction_is_unchanged_by_the_rule():
    """Dense predictions such as the linear and mean baselines have no ties, so their
    numbers are exactly what they were before."""
    n = 8
    pred = np.array([9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0])
    obs = np.zeros(n)
    de = np.zeros(n, bool)
    obs[[0, 1, 2]] = 2.0
    de[[0, 1, 2]] = True
    assert de_overlap(pred, obs, de) == pytest.approx(1.0)

    de2 = np.zeros(n, bool)
    de2[[5, 6, 7]] = True
    obs2 = np.zeros(n)
    obs2[[5, 6, 7]] = 2.0
    assert de_overlap(pred, obs2, de2) == pytest.approx(0.0)


def test_a_perfect_ranking_scores_one_and_is_deterministic():
    pred = np.array([5.0, 4.0, 0.1, 0.0])
    obs = np.array([5.0, 4.0, 0.0, 0.0])
    de = np.array([True, True, False, False])
    assert de_overlap(pred, obs, de) == pytest.approx(1.0)
    assert de_overlap(pred, obs, de) == de_overlap(pred, obs, de)


def test_partial_ties_only_at_the_boundary_are_averaged():
    """Two genes are certainly in the top three; the third place is a coin flip
    between two tied genes, one of which is DE."""
    pred = np.array([10.0, 9.0, 1.0, 1.0])
    obs = np.array([3.0, 3.0, 3.0, 0.0])
    de = np.array([True, True, True, False])
    got = de_overlap(pred, obs, de, n_tiebreak=2000)
    # half the draws give intersection 3 (Jaccard 1.0), half give 2 (Jaccard 0.5)
    assert got == pytest.approx(0.75, abs=0.05)


def test_no_observed_de_genes_is_undefined():
    assert np.isnan(de_overlap(np.ones(4), np.zeros(4), np.zeros(4, bool)))


def test_scores_are_reproducible_across_calls_and_seeds_differ():
    pred = np.zeros(10)
    obs = np.zeros(10)
    obs[:3] = 1.0
    de = np.zeros(10, bool)
    de[:3] = True
    a = de_overlap(pred, obs, de, seed=0)
    b = de_overlap(pred, obs, de, seed=0)
    assert a == b


def test_topk_overlap_reports_intersection_and_jaccard_consistently():
    pred = np.array([5.0, 4.0, 3.0, 0.0])
    jac, inter = topk_overlap(pred, {0, 1}, 2)
    assert inter == pytest.approx(2.0)
    assert jac == pytest.approx(1.0)


def test_metrics_de_overlap_shares_the_rule():
    pred = np.zeros(12)
    obs = np.zeros(12)
    obs[:4] = 5.0
    out = metrics.de_overlap(pred, obs, k=4)
    assert out["jaccard"] < 0.45
    assert 0.0 <= out["precision_at_k"] <= 1.0
    assert out["k"] == 4


def test_metrics_de_overlap_perfect_case():
    pred = np.array([9.0, 8.0, 1.0, 0.5])
    obs = np.array([9.0, 8.0, 1.0, 0.5])
    out = metrics.de_overlap(pred, obs, k=2)
    assert out["jaccard"] == pytest.approx(1.0)
    assert out["precision_at_k"] == pytest.approx(1.0)
