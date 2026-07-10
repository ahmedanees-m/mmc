"""Precondition edge-classification logic, tested offline with synthetic effects
(no store access). Each tuple is (effect_size, fdr, crossguide_r, n_downstream);
the activity n_downstream is the power gate."""
from mmc.data.precondition import ACTIVE_MIN, _label

HI = ACTIVE_MIN + 500   # clearly active
LO = ACTIVE_MIN - 1     # inactive


def edge(effect, fdr, activity):
    return (effect, fdr, None, activity)


def test_conserved_same_sign_significant_both():
    assert _label(edge(-2.0, 1e-4, HI), edge(-1.5, 1e-3, HI)) == "conserved"


def test_rewired_opposite_sign():
    assert _label(edge(-2.0, 1e-4, HI), edge(+2.0, 1e-4, HI)) == "rewired"


def test_rewired_significant_in_one_state_only():
    assert _label(edge(-2.0, 1e-4, HI), edge(-0.1, 0.9, HI)) == "rewired"


def test_no_effect_when_active_but_never_significant():
    assert _label(edge(0.0, 0.9, HI), edge(0.1, 0.8, HI)) == "no_effect"


def test_untestable_when_a_perturbation_is_inactive():
    assert _label(edge(-2.0, 1e-4, LO), edge(-2.0, 1e-4, HI)) == "untestable"


def test_untestable_when_not_measured():
    assert _label(None, edge(-2.0, 1e-4, HI)) == "untestable"
