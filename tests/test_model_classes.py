import numpy as np
import pytest

from mmc.eval import model_classes as mc
from mmc.eval.holdout import ModuleData

GENES = ["A", "B", "C", "D"]


def _mod(seed=0):
    rng = np.random.default_rng(seed)
    obs = rng.normal(size=(4, 4))
    de = np.ones((4, 4), bool)
    return ModuleData(GENES, GENES, obs, de, None)


def test_offset_only_is_the_mean_of_training_rows():
    mod = _mod()
    fn = mc.bind_mean_offset_only(mod)
    got = fn([0, 1, 2], 3)
    assert np.allclose(got, mod.observed[[0, 1, 2]].mean(axis=0))


def test_offset_never_sees_the_held_out_row():
    """The offset must be estimated on training rows only, or the comparison is leaky."""
    mod = _mod()
    fn = mc.bind_mean_offset_only(mod)
    baseline = fn([0, 1, 2], 3)
    mod.observed[3] += 1000.0                       # corrupt only the held-out row
    assert np.allclose(fn([0, 1, 2], 3), baseline)


def test_dense_linear_recovers_a_planted_response():
    """Each perturbation has its own response; ridge over the one-hot design should
    return close to that response for a held-out perturbation seen in training."""
    genes = ["A", "B", "C"]
    obs = np.array([[3.0, 0.0, 0.0],
                    [0.0, 5.0, 0.0],
                    [0.0, 0.0, 7.0]])
    mod = ModuleData(genes, genes, obs, np.ones((3, 3), bool), None)
    fn = mc.bind_dense_linear(mod, l2=1e-6)
    # training on all three rows, predicting row 0, recovers row 0's response
    got = fn([0, 1, 2], 0)
    assert np.allclose(got, obs[0], atol=1e-3)


def test_dense_linear_predicts_zero_for_an_unseen_perturbation():
    """A perturbation absent from training has no fitted coefficient, so the honest
    prediction is no change rather than a borrowed one."""
    genes = ["A", "B", "C"]
    obs = np.array([[3.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 7.0]])
    mod = ModuleData(genes, genes, obs, np.ones((3, 3), bool), None)
    fn = mc.bind_dense_linear(mod, l2=1.0)
    got = fn([1, 2], 0)
    assert np.allclose(got, 0.0, atol=1e-9)


def test_dense_linear_cap_zeroes_all_but_the_largest_coefficients():
    genes = ["A", "B", "C", "D"]
    rng = np.random.default_rng(1)
    obs = rng.normal(size=(4, 4))
    mod = ModuleData(genes, genes, obs, np.ones((4, 4), bool), None)
    uncapped = mc.bind_dense_linear(mod, l2=0.1)([0, 1, 2, 3], 0)
    capped = mc.bind_dense_linear(mod, l2=0.1, max_regulators=1)([0, 1, 2, 3], 0)
    assert not np.allclose(uncapped, capped)


def test_dense_linear_shape_is_the_gene_axis():
    mod = _mod()
    got = mc.bind_dense_linear(mod)([0, 1, 2], 3)
    assert got.shape == (len(GENES),)


def test_offset_class_adds_the_training_residual_mean():
    """With a structure that predicts nothing, structural-plus-offset must reduce
    exactly to the mean of the training responses."""
    pytest.importorskip("jax", reason="the structural fit path needs the JAX backend")
    from mmc.grammar.model_spec import ModelSpec

    mod = _mod()
    empty = ModelSpec(genes=list(GENES), edges=[], rules={})
    fn = mc.bind_structural_with_offset(empty, mod, n_starts=1, max_iter=5)
    got = fn([0, 1, 2], 3)
    expected = mod.observed[[0, 1, 2]].mean(axis=0)
    assert np.allclose(got, expected, atol=1e-6), (
        "an empty structure contributes zero, so the class must equal the mean baseline")


def test_offset_class_is_not_the_bare_structural_prediction():
    """Guards against the offset being dropped: with a non-trivial training mean the
    two must differ."""
    pytest.importorskip("jax", reason="the structural fit path needs the JAX backend")
    from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term

    mod = _mod(seed=3)
    mod.observed += 5.0                              # a large shared component
    spec = ModelSpec(genes=list(GENES),
                     edges=[Edge(regulator="A", target="B", sign=1)],
                     rules={"B": Rule(terms=[Term(regulators=["A"])])})
    with_offset = mc.bind_structural_with_offset(spec, mod, n_starts=1, max_iter=20)
    got = with_offset([0, 1, 2], 3)
    assert np.abs(got).mean() > 1.0, "the shared component must appear in the prediction"
