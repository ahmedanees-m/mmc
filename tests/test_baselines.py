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

