"""The pre-registered splits and the leakage audit."""
from mmc.data.splits import build, leakage_audit


def test_leakage_audit_passes():
    audit = leakage_audit(build())
    assert audit["passed"], audit["checks"]


def test_tierA_is_strict_transfer():
    s = build()
    # Tier A sees no test-state perturbation at fit time.
    assert all(state == s.train_state for state, _ in s.visible("A"))
    assert s.visible("A").isdisjoint(s.scored("A"))


def test_tierB_discovery_and_heldout_are_disjoint():
    s = build()
    assert set(s.tierB_discovery).isdisjoint(s.tierB_heldout)
    # what Tier B sees never overlaps what Tier B scores
    assert s.visible("B").isdisjoint(s.scored("B"))


def test_train_and_incontext_cover_the_regulators_once():
    s = build()
    assert set(s.train_perts).isdisjoint(s.incontext_heldout)
    assert (set(s.train_perts) | set(s.incontext_heldout)) == set(s.tierA_perts)
