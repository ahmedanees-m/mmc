"""Stage-0 module-selection precondition test at the edge level.

For a candidate module (candidate regulators and target genes) and an ordered
condition pair (train, test), classify each regulator-to-target edge as conserved,
rewired, no-effect, or untestable from the per-condition knockdown effects, and
compute the conservation and rewiring fractions. A module and direction pass when
conservation is at least 0.5, rewiring is above 0, and there are at least N_MIN real
(conserved plus rewired) edges. The per-edge classification is the ground-truth
scaffold that Step 6's predicted conserved and rewired map is scored against.

Classification for one regulator-to-target edge, given the perturbation's effect,
FDR, and downstream activity in each condition. Activity, the count of downstream
genes the knockdown moved, is the power gate: it is available for every perturbation
(cross-guide concordance is null for single-guide perturbations) and confirms the
assay worked, so a null under an active perturbation is a trustworthy null.
    active       the perturbation moved at least ACTIVE_MIN downstream genes
    significant  FDR below FDR_SIG (an effect on this target is present)
    conserved    active in both, significant in both, same sign
    rewired      active in both, and either significant with opposite sign, or
                 significant in exactly one (a trustworthy null in the other)
    no_effect    active in both, significant in neither (no edge in either state)
    untestable   inactive in at least one condition, or not measured
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from ..shared import store

FDR_SIG = 0.10
ACTIVE_MIN = 20        # a knockdown that moved at least this many genes is well-powered
CROSSGUIDE_MIN = 0.30  # reported alongside; null for single-guide perturbations
N_MIN = 8

MODULES: dict[str, dict[str, list[str]]] = {
    "Th2_GATA3": {
        "regulators": ["GATA3", "STAT6", "TBX21", "STAT4"],
        "targets": ["GATA3", "STAT6", "TBX21", "STAT4", "IL4", "IL5", "IL13"],
    },
    "TCR_signalosome": {
        "regulators": ["CD3E", "ZAP70", "LAT", "LCP2", "PLCG1", "PRKCQ"],
        "targets": ["ZAP70", "LAT", "LCP2", "PLCG1", "PRKCQ",
                    "IL2", "NFKB1", "RELA", "FOS", "JUN"],
    },
}

DIRECTIONS: list[tuple[str, str]] = [
    ("Rest", "Stim8hr"), ("Rest", "Stim48hr"), ("Stim8hr", "Stim48hr"),
]


@dataclass
class EdgeClass:
    regulator: str
    target: str
    train: tuple | None   # (effect_size, fdr, crossguide_r) or None if not measured
    test: tuple | None
    label: str


def _num(x) -> float | None:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    return float(x)


def _active(n_downstream) -> bool:
    n = _num(n_downstream)
    return n is not None and n >= ACTIVE_MIN


def _significant(fdr) -> bool:
    f = _num(fdr)
    return f is not None and f < FDR_SIG


def _index(df) -> dict[tuple[str, str], tuple]:
    out: dict[tuple[str, str], tuple] = {}
    for _, row in df.iterrows():
        out[(row["perturbation"], row["target_gene"])] = (
            row["effect_size"], row["fdr"], row["crossguide_r"], row["n_downstream"])
    return out


def _label(a: tuple | None, b: tuple | None) -> str:
    if a is None or b is None:
        return "untestable"
    (e_a, f_a, _, nd_a), (e_b, f_b, _, nd_b) = a, b
    if not (_active(nd_a) and _active(nd_b)):
        return "untestable"
    sig_a, sig_b = _significant(f_a), _significant(f_b)
    if sig_a and sig_b:
        return "conserved" if (e_a > 0) == (e_b > 0) else "rewired"
    if sig_a or sig_b:
        return "rewired"
    return "no_effect"


def classify_module(name: str, train_cond: str, test_cond: str) -> dict:
    m = MODULES[name]
    regs, tgts = m["regulators"], m["targets"]
    tr = _index(store.module_effects(regs, tgts, train_cond))
    te = _index(store.module_effects(regs, tgts, test_cond))
    edges: list[EdgeClass] = []
    for reg in regs:
        for tgt in tgts:
            if reg == tgt:
                continue
            a, b = tr.get((reg, tgt)), te.get((reg, tgt))
            edges.append(EdgeClass(reg, tgt, a, b, _label(a, b)))
    counts = {k: sum(1 for e in edges if e.label == k)
              for k in ("conserved", "rewired", "no_effect", "untestable")}
    real = counts["conserved"] + counts["rewired"]
    conservation = counts["conserved"] / real if real else None
    rewiring = counts["rewired"] / real if real else None
    passed = (conservation is not None and conservation >= 0.5
              and counts["rewired"] > 0 and real >= N_MIN)
    return {
        "module": name, "train": train_cond, "test": test_cond,
        "counts": counts, "real_edges": real,
        "conservation": conservation, "rewiring": rewiring, "passed": passed,
        "edges": edges,
    }


def coverage(name: str) -> dict:
    """Which module genes are present in the measured transcriptome."""
    m = MODULES[name]
    genes = sorted(set(m["regulators"]) | set(m["targets"]))
    measured = store.measured_genes()
    return {"present": [g for g in genes if g in measured],
            "missing": [g for g in genes if g not in measured]}


def run_all() -> list[dict]:
    return [classify_module(name, train, test)
            for name in MODULES for train, test in DIRECTIONS]
