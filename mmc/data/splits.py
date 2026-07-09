"""Build the pre-registered splits and enforce the leakage boundaries by construction.

The split definition is fixed in prereg/PREREG.md section 4 and is reproduced here as
data, not re-derived. Training is in the early state (Stim8hr); the decisive test is
transfer to the late state (Stim48hr). Tier A is strict transfer: the frozen structure
predicts every late-state perturbation having seen zero late-state perturbations at fit
time. Tier B is few-shot rewiring discovery: the loop may see a discovery subset of
late-state perturbations, and is scored on a disjoint held-out subset it never sees.

leakage_audit checks these boundaries hold and is a committed deliverable: the score is
only meaningful if the scored perturbations were genuinely unseen when the model was
frozen.
"""
from __future__ import annotations

from dataclasses import dataclass

TRAIN_STATE = "Stim8hr"
TEST_STATE = "Stim48hr"

# PREREG section 4, primary module (TCR signalosome). Structure is trained on four of
# the six regulator knockdowns in the early state; the other two are the in-context
# held-out set for the Step-5 gate. The late state supplies the transfer tests.
_TRAIN_PERTS = ("CD3E", "ZAP70", "LAT", "PRKCQ")
_INCONTEXT_HELDOUT = ("LCP2", "PLCG1")
_TIER_A = ("CD3E", "ZAP70", "LAT", "LCP2", "PLCG1", "PRKCQ")
_TIER_B_DISCOVERY = ("CD3E", "ZAP70", "LAT")
_TIER_B_HELDOUT = ("LCP2", "PLCG1", "PRKCQ")


@dataclass(frozen=True)
class Splits:
    module: str
    train_state: str
    test_state: str
    train_perts: tuple[str, ...]        # structure-training perturbations, train state
    incontext_heldout: tuple[str, ...]  # held out in the train state for the Step-5 gate
    tierA_perts: tuple[str, ...]        # strict transfer: every test-state perturbation
    tierB_discovery: tuple[str, ...]    # test-state perturbations the loop may see
    tierB_heldout: tuple[str, ...]      # test-state perturbations scored, never seen

    def visible(self, tier: str) -> set[tuple[str, str]]:
        """(state, perturbation) pairs the loop and fitter may see for a given tier.

        Tier A sees only train-state perturbations (strict transfer). Tier B additionally
        sees the discovery subset in the test state (few-shot rewiring discovery).
        """
        base = {(self.train_state, p) for p in self.train_perts}
        if tier == "A":
            return base
        if tier == "B":
            return base | {(self.test_state, p) for p in self.tierB_discovery}
        raise ValueError(f"unknown tier {tier!r}")

    def scored(self, tier: str) -> set[tuple[str, str]]:
        """(state, perturbation) pairs that are scored for a given tier."""
        if tier == "A":
            return {(self.test_state, p) for p in self.tierA_perts}
        if tier == "B":
            return {(self.test_state, p) for p in self.tierB_heldout}
        raise ValueError(f"unknown tier {tier!r}")


def build(module: str = "TCR_signalosome") -> Splits:
    if module == "TCR_signalosome":
        return Splits(
            module=module,
            train_state=TRAIN_STATE,
            test_state=TEST_STATE,
            train_perts=_TRAIN_PERTS,
            incontext_heldout=_INCONTEXT_HELDOUT,
            tierA_perts=_TIER_A,
            tierB_discovery=_TIER_B_DISCOVERY,
            tierB_heldout=_TIER_B_HELDOUT,
        )
    return _generic(module)


def _generic(module: str) -> Splits:
    """Derive a leakage-safe split for a module from the perturbations measured in both
    states. Deterministic by sorted gene order: the first 70 percent train structure and
    the rest are the in-context held-out set; the first and second halves are the Tier B
    discovery and held-out subsets. Added as a dated PREREG amendment for modules beyond
    the pre-registered TCR signalosome."""
    from . import module_extract

    tr = set(module_extract.observed_deltas(module, TRAIN_STATE))
    te = set(module_extract.observed_deltas(module, TEST_STATE))
    regs = sorted(tr & te)
    n = len(regs)
    n_train = max(1, round(0.7 * n))
    half = max(1, n // 2)
    return Splits(
        module=module, train_state=TRAIN_STATE, test_state=TEST_STATE,
        train_perts=tuple(regs[:n_train]), incontext_heldout=tuple(regs[n_train:]),
        tierA_perts=tuple(regs), tierB_discovery=tuple(regs[:half]),
        tierB_heldout=tuple(regs[half:]),
    )


def leakage_audit(s: Splits) -> dict:
    """Verify the pre-registered leakage boundaries hold by construction."""
    checks = {
        # Tier A is strict: nothing in the test state is visible at fit time.
        "tierA_sees_no_test_state": all(st != s.test_state for st, _ in s.visible("A")),
        # what is scored is never in what was visible, per tier.
        "tierA_scored_unseen": s.visible("A").isdisjoint(s.scored("A")),
        "tierB_scored_unseen": s.visible("B").isdisjoint(s.scored("B")),
        # the Tier B discovery and held-out subsets are disjoint.
        "tierB_discovery_heldout_disjoint":
            set(s.tierB_discovery).isdisjoint(s.tierB_heldout),
        # the in-context gate holds out perturbations the structure was not trained on.
        "train_incontext_disjoint": set(s.train_perts).isdisjoint(s.incontext_heldout),
        # the train and in-context sets together cover the regulator set exactly once.
        "train_plus_incontext_cover_regulators":
            (set(s.train_perts) | set(s.incontext_heldout)) == set(s.tierA_perts),
    }
    return {"passed": all(checks.values()), "checks": checks}
