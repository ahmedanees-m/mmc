import numpy as np
import pytest

from mmc.eval import compare
from mmc.eval.holdout import ModuleData


def _mod(n_perts=6, n_genes=5, seed=0):
    rng = np.random.default_rng(seed)
    obs = rng.normal(size=(n_perts, n_genes))
    de = np.zeros((n_perts, n_genes), bool)
    de[:, :2] = True
    return ModuleData([f"G{j}" for j in range(n_genes)],
                      [f"P{i}" for i in range(n_perts)], obs, de, None)


def test_paired_delta_recovers_a_constant_offset():
    a = np.array([0.5, 0.6, 0.7, 0.8])
    b = a - 0.1
    r = compare.paired_delta(a, b, n_boot=2000, seed=0)
    assert r["delta"] == pytest.approx(0.1, abs=1e-9)
    # a constant within-fold offset has zero paired variance, so the interval collapses
    assert r["lo"] == pytest.approx(0.1, abs=1e-9)
    assert r["advantage"] is True
    assert r["n_folds"] == 4


def test_paired_delta_is_not_fooled_by_shared_between_fold_variance():
    """The reason the paired form replaces the overlapping-interval check.

    Both arms swing hugely across folds and one is uniformly better. Marginal
    intervals overlap; the paired interval does not.
    """
    b = np.array([0.05, 0.95, 0.10, 0.90, 0.15, 0.85])
    a = b + 0.04
    assert compare.paired_delta(a, b, n_boot=4000, seed=0)["advantage"] is True
    marg_a = compare.bootstrap_mean(a, n_boot=4000, seed=0)
    marg_b = compare.bootstrap_mean(b, n_boot=4000, seed=0)
    assert marg_a["lo"] < marg_b["hi"]      # the marginal check would miss it


def test_paired_delta_no_advantage_when_interval_spans_zero():
    rng = np.random.default_rng(1)
    a = rng.normal(0.5, 0.2, size=12)
    b = rng.normal(0.5, 0.2, size=12)
    r = compare.paired_delta(a, b, n_boot=4000, seed=0)
    assert r["lo"] <= 0 <= r["hi"]
    assert r["advantage"] is False


def test_paired_delta_drops_undefined_folds_from_both_arms():
    a = np.array([0.4, np.nan, 0.6, 0.8])
    b = np.array([0.1, 0.2, np.nan, 0.5])
    r = compare.paired_delta(a, b, n_boot=500, seed=0)
    assert r["n_folds"] == 2                       # folds 0 and 3 only
    assert r["delta"] == pytest.approx(((0.4 - 0.1) + (0.8 - 0.5)) / 2)


def test_paired_delta_rejects_misaligned_arms():
    with pytest.raises(ValueError):
        compare.paired_delta(np.zeros(3), np.zeros(4))


def test_paired_delta_all_nan_is_reported_not_raised():
    r = compare.paired_delta(np.array([np.nan]), np.array([np.nan]))
    assert r["n_folds"] == 0 and r["advantage"] is False


def test_benjamini_hochberg_matches_the_published_worked_example():
    """Benjamini and Hochberg 1995, Table 1: at q = 0.05 the first four reject."""
    p = [0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344,
         0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.000]
    assert compare.benjamini_hochberg(p, q=0.05) == [True] * 4 + [False] * 11


def test_benjamini_hochberg_is_step_up_not_per_test():
    """0.04 fails its own threshold but is rescued because a larger rank passes."""
    p = [0.001, 0.04, 0.045]
    assert compare.benjamini_hochberg(p, q=0.05) == [True, True, True]
    # per-test comparison would have rejected only the first
    assert [x <= 0.05 * (i + 1) / 3 for i, x in enumerate(sorted(p))] == [True, False, True]


def test_benjamini_hochberg_preserves_input_order():
    assert compare.benjamini_hochberg([0.9, 0.0001], q=0.05) == [False, True]
    assert compare.benjamini_hochberg([]) == []


def test_permutation_p_is_small_for_a_consistent_gap():
    a = np.full(10, 0.6)
    b = np.full(10, 0.4)
    assert compare.permutation_p(a, b, n_perm=2000, seed=0) < 0.01


def test_permutation_p_is_uninformative_for_identical_arms():
    a = np.linspace(0.1, 0.9, 10)
    assert compare.permutation_p(a, a.copy(), n_perm=1000, seed=0) > 0.2


def test_baseline_sources_use_only_training_rows():
    mod = _mod()
    n = len(mod.perts)
    mean_fn = compare.bind_mean(mod)
    for i in range(n):
        train = [j for j in range(n) if j != i]
        got = mean_fn(train, i)
        assert np.allclose(got, mod.observed[train].mean(axis=0))
    assert np.allclose(compare.bind_zero(mod)([0, 1], 2), 0.0)


def test_fold_predictions_shape_and_scoring():
    mod = _mod()
    preds = compare.fold_predictions(mod, compare.bind_mean(mod))
    assert preds.shape == mod.observed.shape
    scores = compare.score_predictions(mod, preds)
    assert set(scores) == {"de_overlap", "acc_deg"}
    assert scores["de_overlap"].shape == (len(mod.perts),)


def test_comparator_table_orders_by_metric_and_flags_advantage():
    mod = _mod()
    results = {}
    for name in ("linear", "better", "worse"):
        preds = compare.fold_predictions(mod, compare.bind_mean(mod))
        results[name] = compare.SourceResult(name, preds, compare.score_predictions(mod, preds))
    # force a known ordering on the scores
    results["better"].scores["de_overlap"] = results["linear"].scores["de_overlap"] + 0.2
    results["worse"].scores["de_overlap"] = results["linear"].scores["de_overlap"] - 0.2
    rows = compare.comparator_table(results, reference="linear")
    assert [r["source"] for r in rows] == ["better", "linear", "worse"]
    assert rows[0]["advantage"] is True
    assert rows[2]["advantage"] is False


def test_jaccard_edges():
    from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term

    def spec(pairs):
        rules = {}
        for r, t in pairs:
            rules.setdefault(t, []).append(r)
        return ModelSpec(genes=["A", "B", "C"],
                         edges=[Edge(regulator=r, target=t, sign=1) for r, t in pairs],
                         rules={t: Rule(terms=[Term(regulators=rs)]) for t, rs in rules.items()})

    a = spec([("A", "B"), ("A", "C")])
    b = spec([("A", "B")])
    assert compare.jaccard_edges(a, a) == 1.0
    assert compare.jaccard_edges(a, b) == 0.5
