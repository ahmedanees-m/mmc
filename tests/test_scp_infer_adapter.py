import numpy as np
import pytest

from mmc.compare import scp_infer_adapter as s4
from mmc.eval.holdout import ModuleData
from mmc.grammar.model_spec import MAX_REGS_PER_TERM

sklearn = pytest.importorskip("sklearn")

GENES = ["A", "B", "C", "D", "E"]
PERTS = ["A", "B", "C"]


def _mod(seed=0):
    """A module where knocking down A lowers D and raises E, and B lowers C."""
    rng = np.random.default_rng(seed)
    obs = 0.01 * rng.normal(size=(3, 5))
    de = np.zeros((3, 5), bool)
    gi = {g: i for i, g in enumerate(GENES)}
    obs[0, gi["D"]] = -3.0      # A knockdown lowers D, so A activates D
    obs[0, gi["E"]] = +2.5      # A knockdown raises E, so A represses E
    obs[1, gi["C"]] = -2.0      # B activates C
    de[0, gi["D"]] = de[0, gi["E"]] = de[1, gi["C"]] = True
    return ModuleData(GENES, PERTS, obs, de, None)


def test_mean_difference_recovers_signs_from_the_knockdown_direction():
    edges = s4.mean_difference_edges(_mod())
    assert ("A", "D", 1) in edges, "a target that falls on knockdown is activated"
    assert ("A", "E", -1) in edges, "a target that rises on knockdown is repressed"
    assert ("B", "C", 1) in edges


def test_mean_difference_ranks_by_effect_size():
    edges = s4.mean_difference_edges(_mod())
    assert edges[0][:2] == ("A", "D")          # |-3.0| is the largest effect


def test_mean_difference_uses_only_de_entries_by_default():
    edges = s4.mean_difference_edges(_mod())
    assert len(edges) == 3
    loose = s4.mean_difference_edges(_mod(), require_de=False)
    assert len(loose) > len(edges)


def test_mean_difference_never_emits_a_self_edge():
    edges = s4.mean_difference_edges(_mod(), require_de=False, max_edges=99)
    assert all(r != t for r, t, _ in edges)


def test_edge_sets_are_grammar_legal():
    for fn in (s4.mean_difference_edges, s4.grnboost2_edges):
        edges = fn(_mod(), max_edges=99)
        in_deg: dict[str, int] = {}
        pairs = set()
        for r, t, sign in edges:
            assert sign in (1, -1)
            assert (r, t) not in pairs
            pairs.add((r, t))
            in_deg[t] = in_deg.get(t, 0) + 1
        assert not in_deg or max(in_deg.values()) <= MAX_REGS_PER_TERM


def test_max_edges_is_respected():
    assert len(s4.mean_difference_edges(_mod(), require_de=False, max_edges=2)) == 2


def test_grnboost2_returns_signed_edges_from_perturbed_regulators():
    edges = s4.grnboost2_edges(_mod(), max_edges=10)
    assert edges
    assert all(r in PERTS for r, _t, _s in edges)
    assert all(s in (1, -1) for _r, _t, s in edges)


def test_grnboost2_is_deterministic_under_a_seed():
    a = s4.grnboost2_edges(_mod(), seed=3, n_estimators=50)
    b = s4.grnboost2_edges(_mod(), seed=3, n_estimators=50)
    assert a == b


def test_grnboost2_handles_a_module_with_too_few_regulators():
    mod = ModuleData(["A", "B"], ["A"], np.zeros((1, 2)), np.zeros((1, 2), bool), None)
    assert s4.grnboost2_edges(mod) == []


def test_unavailable_methods_explain_why_rather_than_failing_silently():
    for method in ("gies", "dcdi-g", "bicycle", "avici", "sdcd"):
        with pytest.raises(ValueError, match="per-cell interventional data"):
            s4.algorithmic_edges(_mod(), method)


def test_dispatch_reaches_the_runnable_methods():
    for method in s4.METHODS:
        assert isinstance(s4.algorithmic_edges(_mod(), method), list)
