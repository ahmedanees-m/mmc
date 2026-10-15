"""Multi-start fitting and the fit-vs-structure gate on small synthetic circuits.
Skipped where cma is not installed."""
import numpy as np
import pytest

pytest.importorskip("cma")

from mmc.compile.perturb import knockdown
from mmc.compile.simulate import steady_state
from mmc.fit.diagnose import diagnose
from mmc.fit.fit_params import _gene_index, multi_fit
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term


def _chain():
    """A activates B activates C, with a set of true parameters."""
    spec = ModelSpec(
        genes=["A", "B", "C"],
        edges=[Edge(regulator="A", target="B", sign=1),
               Edge(regulator="B", target="C", sign=1)],
        rules={"B": Rule(terms=[Term(regulators=["A"])]),
               "C": Rule(terms=[Term(regulators=["B"])])},
    )
    params = {
        "basal": np.array([1.0, 0.05, 0.05]),
        "decay": np.array([1.0, 1.0, 1.0]),
        "terms": {"B": [{"prod": 3.0, "w": {"A": 3.0}, "theta": {"A": 1.0}}],
                  "C": [{"prod": 3.0, "w": {"B": 3.0}, "theta": {"B": 1.0}}]},
    }
    return spec, params


def _observed_from_model(spec, params):
    gi = _gene_index(spec)
    wt = steady_state(spec, params)
    obs = {}
    for reg in spec.genes:
        d = knockdown(spec, params, reg, wt=wt)
        obs[reg] = {g: float(d[gi[g]]) for g in spec.genes if g != reg}
    return obs


def test_fit_recovers_a_consistent_structure():
    spec, params = _chain()
    observed = _observed_from_model(spec, params)
    fits = multi_fit(spec, observed, n_starts=4, max_iter=60)
    assert fits[0]["loss"] < 0.3          # a consistent structure fits well


def test_diagnose_labels_parametric_when_a_fit_exists():
    spec, params = _chain()
    observed = _observed_from_model(spec, params)
    _best, labels, stats = diagnose(spec, observed, n_starts=4, max_iter=60)
    assert labels and all(v == "parametric" for v in labels.values())
    assert stats["loss_best"] < 0.3


def test_diagnose_flags_a_wrong_sign_as_structural():
    # the structure says A activates B, so a knockdown of A must lower B; the data has
    # B rising on the A knockdown, which no parameter set can reproduce
    spec = ModelSpec(
        genes=["A", "B"],
        edges=[Edge(regulator="A", target="B", sign=1)],
        rules={"B": Rule(terms=[Term(regulators=["A"])])},
    )
    observed = {"A": {"B": +2.0}}
    _best, labels, _stats = diagnose(spec, observed, n_starts=4, max_iter=50)
    assert labels[("A", "B")] == "structural"
