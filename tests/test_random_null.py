import numpy as np
import pytest

from mmc.eval import random_null
from mmc.grammar.model_spec import MAX_REGS_PER_TERM

GENES = [f"G{i}" for i in range(10)]
PERTS = GENES[:6]


def test_sample_spec_respects_the_grammar_bounds():
    rng = np.random.default_rng(0)
    for _ in range(20):
        spec = random_null.sample_spec(GENES, PERTS, 12, rng)
        assert len(spec.edges) == 12
        in_deg: dict[str, int] = {}
        for e in spec.edges:
            assert e.regulator != e.target
            assert e.regulator in PERTS
            assert e.sign in (1, -1)
            in_deg[e.target] = in_deg.get(e.target, 0) + 1
        assert max(in_deg.values()) <= MAX_REGS_PER_TERM
        # no duplicate regulator -> target pair
        pairs = [(e.regulator, e.target) for e in spec.edges]
        assert len(pairs) == len(set(pairs))


def test_sample_spec_is_deterministic_under_a_seed():
    a = random_null.sample_spec(GENES, PERTS, 10, np.random.default_rng(7))
    b = random_null.sample_spec(GENES, PERTS, 10, np.random.default_rng(7))
    assert [(e.regulator, e.target, e.sign) for e in a.edges] == \
           [(e.regulator, e.target, e.sign) for e in b.edges]


def test_sample_spec_varies_across_seeds():
    a = random_null.sample_spec(GENES, PERTS, 10, np.random.default_rng(1))
    b = random_null.sample_spec(GENES, PERTS, 10, np.random.default_rng(2))
    assert {(e.regulator, e.target) for e in a.edges} != \
           {(e.regulator, e.target) for e in b.edges}


def test_sample_spec_fails_loudly_when_the_grammar_cannot_hold_the_request():
    rng = np.random.default_rng(0)
    # 3 targets x in-degree 3 = 9 possible edges at most
    with pytest.raises(ValueError, match="grammar cannot hold"):
        random_null.sample_spec(["A", "B", "C"], ["A", "B"], 50, rng)


def test_sample_spec_rejects_empty_inputs():
    with pytest.raises(ValueError):
        random_null.sample_spec([], [], 1, np.random.default_rng(0))


def test_gated_sampler_keeps_the_edge_set_and_may_split_terms():
    rng = np.random.default_rng(3)
    spec = random_null.sample_spec_gated(GENES, PERTS, 15, rng, p_and=1.0)
    assert len(spec.edges) == 15
    # with p_and = 1 every multi-regulator target should have been split
    multi = [t for t, r in spec.rules.items() if sum(len(x.regulators) for x in r.terms) >= 2]
    assert any(len(spec.rules[t].terms) == 2 for t in multi)


def test_percentile_of():
    null = [0.1, 0.2, 0.3, 0.4]
    assert random_null.percentile_of(0.35, null) == 75.0
    assert random_null.percentile_of(0.0, null) == 0.0
    assert random_null.percentile_of(1.0, null) == 100.0
    assert np.isnan(random_null.percentile_of(float("nan"), null))
    assert np.isnan(random_null.percentile_of(0.5, [float("nan")]))


def test_summarise_null_reports_the_upper_tail():
    v = list(np.linspace(0, 1, 101))
    s = random_null.summarise_null(v)
    assert s["n"] == 101
    assert s["p50"] == pytest.approx(0.5)
    assert s["p95"] == pytest.approx(0.95, abs=1e-9)
    assert s["max"] == pytest.approx(1.0)
    assert random_null.summarise_null([])["n"] == 0
