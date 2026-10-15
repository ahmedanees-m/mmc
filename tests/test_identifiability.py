import numpy as np
import pytest

from mmc.eval import identifiability as ident
from mmc.eval.holdout import ModuleData


def test_effective_rank_is_one_for_a_single_shared_direction():
    """Every perturbation moves the transcriptome the same way, up to scale."""
    direction = np.array([1.0, -2.0, 0.5, 3.0])
    scales = np.array([1.0, 2.0, 3.0, 4.0, 5.0])[:, None]
    r = scales * direction
    assert ident.effective_rank(r) == pytest.approx(1.0, abs=1e-6)


def test_effective_rank_rises_with_independent_directions():
    r = np.eye(5) * 3.0
    assert ident.effective_rank(r) > 3.0


def test_effective_rank_needs_at_least_two_rows():
    assert np.isnan(ident.effective_rank(np.ones((1, 4))))


def test_effective_rank_of_a_degenerate_matrix_is_nan():
    assert np.isnan(ident.effective_rank(np.zeros((4, 3))))


def test_leading_pc_fraction_is_high_when_the_held_out_row_lies_on_pc1():
    train = np.array([[1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [-1.0, 0.0]])
    held = np.array([5.0, 0.0])
    assert ident.leading_pc_variance_fraction(train, held) == pytest.approx(1.0, abs=1e-6)


def test_leading_pc_fraction_is_low_when_the_held_out_row_is_orthogonal():
    train = np.array([[1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [-1.0, 0.0]])
    held = train.mean(axis=0) + np.array([0.0, 5.0])
    assert ident.leading_pc_variance_fraction(train, held) == pytest.approx(0.0, abs=1e-6)


def test_mean_leading_pc_fraction_runs_over_folds():
    rng = np.random.default_rng(0)
    r = rng.normal(size=(6, 4))
    v = ident.mean_leading_pc_fraction(r)
    assert 0.0 <= v <= 1.0
    assert np.isnan(ident.mean_leading_pc_fraction(r[:2]))


def test_perturbation_specific_ratio_is_zero_when_every_response_is_identical():
    r = np.tile(np.array([1.0, -1.0, 2.0]), (5, 1))
    assert ident.perturbation_specific_ratio(r) == pytest.approx(0.0, abs=1e-12)


def test_perturbation_specific_ratio_grows_with_idiosyncratic_responses():
    shared = np.tile(np.array([1.0, 1.0, 1.0]), (6, 1))
    rng = np.random.default_rng(1)
    weak = shared + 0.05 * rng.normal(size=shared.shape)
    strong = shared + 2.0 * rng.normal(size=shared.shape)
    assert ident.perturbation_specific_ratio(weak) < ident.perturbation_specific_ratio(strong)


def test_equivalence_width_spans_the_near_optimal_class():
    # three structures fit the training data equally well and predict very differently
    train_loss = [1.0, 1.02, 1.04, 5.0]
    holdout = [0.10, 0.55, 0.30, 0.99]
    w = ident.equivalence_width(train_loss, holdout, epsilon=1.05)
    assert w["n_in_class"] == 3          # the 5.0 loss is outside the band
    assert w["width"] == pytest.approx(0.45)
    assert w["holdout_max"] == pytest.approx(0.55)
    assert w["best_train_loss"] == pytest.approx(1.0)


def test_equivalence_width_is_zero_when_the_class_agrees():
    w = ident.equivalence_width([1.0, 1.01], [0.4, 0.4], epsilon=1.05)
    assert w["width"] == pytest.approx(0.0)
    assert w["sd"] == pytest.approx(0.0)


def test_equivalence_width_handles_an_empty_trace():
    assert ident.equivalence_width([], [])["n_in_class"] == 0
    assert ident.equivalence_width([np.nan], [np.nan])["n_in_class"] == 0


def test_sign_stability_flags_an_edge_the_data_cannot_orient():
    stable = ("A", "B", 1)
    sets = [
        {stable, ("C", "D", 1)},
        {stable, ("C", "D", -1)},
        {stable, ("C", "D", 1)},
        {stable, ("C", "D", -1)},
    ]
    out = ident.sign_stability(sets, [1.0, 1.0, 1.0, 1.0])
    assert out["n_edges"] == 2
    assert out["per_edge_agreement"]["A->B"] == 1.0
    assert out["per_edge_agreement"]["C->D"] == pytest.approx(0.5)
    assert out["fraction_stable"] == pytest.approx(0.5)


def test_sign_stability_only_looks_inside_the_near_optimal_band():
    sets = [{("A", "B", 1)}, {("A", "B", -1)}]
    # the second structure fits far worse, so it is outside the band and ignored
    out = ident.sign_stability(sets, [1.0, 10.0])
    assert out["per_edge_agreement"]["A->B"] == 1.0
    assert out["fraction_stable"] == 1.0


def test_sign_stability_handles_empty_input():
    assert ident.sign_stability([], [])["n_edges"] == 0


def test_effect_size_summary_counts_de_entries():
    obs = np.array([[2.0, 0.0, -3.0], [0.5, 1.0, 0.0]])
    de = np.array([[True, False, True], [False, True, False]])
    s = ident.effect_size_summary(obs, de)
    assert s["n_de_entries"] == 3
    assert s["de_per_pert_mean"] == pytest.approx(1.5)
    assert s["de_per_pert_min"] == 1
    assert s["n_perts"] == 2 and s["n_genes"] == 3


def test_diagnostics_bundles_the_data_side_quantities():
    rng = np.random.default_rng(2)
    obs = rng.normal(size=(6, 5))
    de = np.zeros((6, 5), bool)
    de[:, :2] = True
    mod = ModuleData([f"G{i}" for i in range(5)], [f"P{i}" for i in range(6)], obs, de, None)
    d = ident.diagnostics(mod)
    assert {"effective_rank", "leading_pc_fraction",
            "perturbation_specific_ratio", "effect_sizes"} <= set(d)
    assert "equivalence_width" not in d          # no trace supplied


def test_diagnostics_picks_up_the_trace_when_given_one():
    rng = np.random.default_rng(3)
    obs = rng.normal(size=(4, 3))
    de = np.ones((4, 3), bool)
    mod = ModuleData(["A", "B", "C"], ["A", "B", "C", "D"][:4], obs, de, None)
    trace = {"train_loss": [1.0, 1.01], "holdout": [0.2, 0.4],
             "edge_sets_raw": [{("A", "B", 1)}, {("A", "B", 1)}]}
    d = ident.diagnostics(mod, trace=trace)
    assert d["equivalence_width"]["n_in_class"] == 2
    assert d["sign_stability"]["fraction_stable"] == 1.0
