"""The negative and positive controls."""
import numpy as np

from mmc.eval import controls


def test_wt_positive_is_perfect():
    obs = {"A": {"G1": 1.0, "G2": -2.0, "G3": 0.5},
           "B": {"G1": -1.0, "G2": 2.0, "G3": 0.0}}
    assert abs(controls.wt_positive(obs, ["G1", "G2", "G3"]) - 1.0) < 1e-9


def test_shuffled_negative_does_not_stay_perfect():
    genes = [f"G{i}" for i in range(10)]
    rng = np.random.default_rng(0)
    obs = {p: {g: float(rng.normal()) for g in genes} for p in ("A", "B", "C", "D")}
    pred = {p: dict(obs[p]) for p in obs}                 # perfect predictions
    neg = controls.shuffled_negative(pred, obs, genes, seed=2)
    # matched to the wrong perturbation, the perfect predictions no longer correlate
    assert neg is None or abs(neg) < 0.9
