"""The mean and linear baselines."""
import numpy as np

from mmc.baselines import consensus, linear, mean

TRAIN = {"A": {"X": 1.0, "Y": -2.0}, "B": {"X": 3.0, "Y": 0.0}}
GENES = ["X", "Y", "Z"]


def test_mean_profile_and_independent_copies():
    prof = mean.fit(TRAIN, GENES)
    assert prof["X"] == 2.0      # mean of 1 and 3
    assert prof["Y"] == -1.0     # mean of -2 and 0
    assert prof["Z"] == 0.0      # unmeasured gene falls back to zero
    pred = mean.predict(prof, ["H1", "H2"])
    assert pred["H1"] == prof and pred["H1"] is not pred["H2"]


def test_persistence_returns_the_train_state_effect():
    pred = linear.persistence(TRAIN, ["A", "B"])
    assert pred["A"] == {"X": 1.0, "Y": -2.0}


def test_reconstruct_predicts_a_finite_vector_for_a_held_out_perturbation():
    train_deltas = {"A": {"B": 1.0, "C": -1.0}, "D": {"B": 0.5, "C": 2.0}}
    genes = ["A", "B", "C", "D"]
    pred = linear.reconstruct(train_deltas, genes, ["A", "D"], ["B"], l2=1.0)
    assert set(pred["B"]) == set(genes)
    assert all(np.isfinite(v) for v in pred["B"].values())


def test_consensus_leans_toward_the_similar_neighbour():
    # in the training state the held-out gene H behaves like discovery gene D1, not D2
    train = {"S1": {"H": 1.0, "D1": 1.0, "D2": -1.0},
             "S2": {"H": 2.0, "D1": 2.0, "D2": -2.0}}
    test_visible = {"D1": {"G": 5.0}, "D2": {"G": -5.0}}
    pred = consensus.predict(train, test_visible, ["G"], ["H"], ["D1", "D2"], ["S1", "S2"])
    assert pred["H"]["G"] > 0        # weighted toward D1, the neighbour H resembles



def _rank_of_predictions(rank):
    """Fit the reduced-rank map on a rank-2 response matrix and return the predictions."""
    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(12)]
    basis = rng.normal(size=(2, len(genes)))
    perts = genes[:9]           # a perturbation is identified by the gene it targets
    loads = rng.normal(size=(len(perts), 2))
    obs = loads @ basis                                  # exactly rank 2 by construction
    train = {p: {g: float(obs[i, j]) for j, g in enumerate(genes)}
             for i, p in enumerate(perts)}
    pred = linear.reconstruct_reduced_rank(train, genes, perts, perts, rank=rank, l2=1e-6)
    return np.vstack([[pred[p][g] for g in genes] for p in perts])


def test_reduced_rank_truncation_holds_the_requested_rank():
    for k in (1, 2):
        s = np.linalg.svd(_rank_of_predictions(k), compute_uv=False)
        assert np.sum(s > 1e-8 * s[0]) <= k


def test_reduced_rank_at_the_true_rank_matches_the_full_ridge():
    # the responses are exactly rank 2, so truncating at 2 should cost essentially nothing
    full = _rank_of_predictions(12)
    at_true = _rank_of_predictions(2)
    assert np.allclose(full, at_true, atol=1e-6)
    # and truncating below the true rank must cost something, or the test proves nothing
    assert not np.allclose(full, _rank_of_predictions(1), atol=1e-6)
