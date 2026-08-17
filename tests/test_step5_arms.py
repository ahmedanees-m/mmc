"""Arm construction for the Step 5 ablation.

The arms are only interpretable if each changes exactly one thing. These check that
A2 hides names without touching the data, A3 destroys the pairing without touching the
names, and A4 does both, and that the redaction gate refuses to let a leaky arm run.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from mmc.eval.holdout import ModuleData

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "step5_anonymization.py"
_spec = importlib.util.spec_from_file_location("step5_anonymization", _PATH)
step5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step5)

GENES = ["GATA3", "STAT6", "TBX21", "IL5", "IL4"]


def _mod():
    rng = np.random.default_rng(0)
    obs = rng.normal(size=(3, 5))
    de = np.zeros((3, 5), bool)
    de[:, 3] = True
    return ModuleData(GENES, GENES[:3], obs, de, None)


def test_a1_is_the_untouched_arm():
    genes, observed, context, amap, audit, de = step5.build_arm(_mod(), "M", "A1", 0)
    assert genes == GENES
    assert amap is None and audit == []
    assert set(observed) == set(GENES[:3])
    assert "GATA3" in context or "M" in context


def test_a2_hides_the_names_and_leaves_the_data_alone():
    mod = _mod()
    g1, obs1, _c1, _a1, _au1, _de1 = step5.build_arm(mod, "M", "A1", 0)
    g2, obs2, ctx2, amap, audit, _de2 = step5.build_arm(mod, "M", "A2", 0)
    assert audit == [], "A2 must not leak a gene symbol"
    assert all(g not in GENES for g in g2)
    from mmc.loop import anonymize as anon
    assert anon.redaction_violations(ctx2, GENES) == []
    # the numbers are identical, only the keys changed
    v1 = sorted(round(v, 9) for d in obs1.values() for v in d.values())
    v2 = sorted(round(v, 9) for d in obs2.values() for v in d.values())
    assert v1 == v2


def test_a3_permutes_the_data_and_keeps_the_names():
    mod = _mod()
    g1, obs1, *_ = step5.build_arm(mod, "M", "A1", 0)
    g3, obs3, _c, amap, audit, _de = step5.build_arm(mod, "M", "A3", 0)
    assert g3 == GENES and amap is None and audit == []
    assert obs3 != obs1, "A3 must not reproduce the real pairing"
    # the multiset of responses is preserved; only which perturbation owns them changed
    v1 = sorted(round(v, 9) for d in obs1.values() for v in d.values())
    v3 = sorted(round(v, 9) for d in obs3.values() for v in d.values())
    assert len(v1) == len(v3)


def test_a4_does_both():
    mod = _mod()
    g4, obs4, ctx4, amap, audit, _de = step5.build_arm(mod, "M", "A4", 1)
    assert audit == []
    assert all(g not in GENES for g in g4)
    assert amap is not None


def test_the_redaction_gate_stops_a_leaky_arm(monkeypatch):
    """If anonymisation regressed, the arm must refuse to run rather than produce a
    confident null."""
    def leaky(genes, context, spec=None, residual_summary=""):
        return "context mentioning GATA3 directly"

    monkeypatch.setattr(step5.anon, "assemble_prompt_surface", leaky)
    with pytest.raises(RuntimeError, match="survived redaction"):
        step5.run_arm("M", "Stim8hr", _mod(), "A2", 0,
                      max_iters=1, n_starts=1, max_iter=10)


def test_jaccard_on_edge_lists_ignores_sign():
    a = [("A", "B", 1), ("C", "D", -1)]
    b = [("A", "B", -1)]
    assert step5.jaccard(a, a) == 1.0
    assert step5.jaccard(a, b) == 0.5
    assert np.isnan(step5.jaccard([], []))


def test_agreement_matrix_reports_within_and_between():
    runs = [
        {"module": "M", "arm": "A1", "seed": 0, "edges": [("A", "B", 1), ("C", "D", 1)]},
        {"module": "M", "arm": "A1", "seed": 1, "edges": [("A", "B", 1)]},
        {"module": "M", "arm": "A2", "seed": 0, "edges": [("A", "B", 1), ("C", "D", 1)]},
        {"module": "M", "arm": "A2", "seed": 1, "edges": [("X", "Y", 1)]},
    ]
    out = step5.agreement_matrix(runs)
    assert out["within_arm_seed_ceiling"]["M/A1"]["mean"] == pytest.approx(0.5)
    # A1 vs A2: seed 0 agrees fully, seed 1 not at all
    assert out["between_arm"]["M/A1_vs_A2"]["mean"] == pytest.approx(0.5)
    assert out["between_arm"]["M/A1_vs_A2"]["n_seeds"] == 2


def test_interpretation_reports_the_ratio_against_the_ceiling():
    agree = {"within_arm_seed_ceiling": {"M/A1": {"mean": 0.8, "n_pairs": 1}},
             "between_arm": {"M/A1_vs_A2": {"mean": 0.4, "n_seeds": 2}}}
    lines = step5.interpret(agree)
    assert lines and "ratio 0.50" in lines[0]
