"""Mechanical module extraction for the Step 7 scale-up (PREREG_v4 section 8).

Section 8 fixes the selection procedure and requires that **every module extracted is
reported, including those that fail power screening**, so the set cannot be curated
toward the conclusion. That requirement shapes this module: `generate` returns every
candidate with its screen verdict attached, and it is the caller's job to run the ones
that pass, not this code's job to hide the ones that do not.

Two candidate sources, as the plan specifies.

*Annotation-driven.* A transcription factor's annotated regulon, intersected with the
genes the atlas measures and perturbs. This is the regime where a regulatory structure
should have the best chance, because the module is defined by regulation.

*Data-driven.* Genes clustered by the correlation of their response profiles across
perturbations. This is the regime the Step 3 calibration says the data actually lives
in, since the response matrix is low rank; a co-response cluster is a shared program,
not a causal neighbourhood, and the contrast between the two sources is itself
informative for the regime map.

The screen is the existing power precondition, not a new one: enough perturbations,
enough DE entries, and enough folds carrying a DE gene for the paired statistic to have
anything to work with. Th2_GATA3 ran with two scoreable folds and produced intervals
spanning the whole metric, which is what the fold-count floor exists to prevent.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

MIN_PERTS = 8
MIN_DE_ENTRIES = 30
MIN_SCOREABLE_FOLDS = 8
MAX_GENES = 40          # section 2.5 keeps modules under the DCDI scalability limit


@dataclass
class Candidate:
    """A proposed module and the screen's verdict on it."""

    name: str
    source: str                       # "regulon" or "coresponse"
    genes: list[str]
    seed_gene: str = ""
    n_perts: int = 0
    n_de_entries: int = 0
    n_scoreable_folds: int = 0
    passed: bool = False
    reasons: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {"name": self.name, "source": self.source, "seed_gene": self.seed_gene,
                "n_genes": len(self.genes), "genes": list(self.genes),
                "n_perts": self.n_perts, "n_de_entries": self.n_de_entries,
                "n_scoreable_folds": self.n_scoreable_folds,
                "passed": self.passed, "reasons": list(self.reasons)}


def screen(cand: Candidate, *, min_perts: int = MIN_PERTS,
           min_de: int = MIN_DE_ENTRIES,
           min_folds: int = MIN_SCOREABLE_FOLDS) -> Candidate:
    """Apply the power precondition and record why a candidate failed."""
    reasons = []
    if cand.n_perts < min_perts:
        reasons.append(f"only {cand.n_perts} perturbations, need {min_perts}")
    if cand.n_de_entries < min_de:
        reasons.append(f"only {cand.n_de_entries} DE entries, need {min_de}")
    if cand.n_scoreable_folds < min_folds:
        reasons.append(f"only {cand.n_scoreable_folds} scoreable folds, need {min_folds}")
    if len(cand.genes) > MAX_GENES:
        reasons.append(f"{len(cand.genes)} genes exceeds the {MAX_GENES} cap")
    cand.reasons = reasons
    cand.passed = not reasons
    return cand


def characterise(genes: list[str], observed: np.ndarray, de_mask: np.ndarray) -> dict:
    """Shape figures the screen needs, from an already-built response matrix."""
    per_pert = de_mask.sum(axis=1)
    return {"n_perts": int(observed.shape[0]),
            "n_de_entries": int(de_mask.sum()),
            "n_scoreable_folds": int((per_pert > 0).sum())}


def regulon_candidates(regulon, measured: set[str], perturbed: set[str], *,
                       min_targets: int = 6, max_genes: int = MAX_GENES,
                       max_modules: int = 40) -> list[Candidate]:
    """One candidate per transcription factor with enough measured, perturbed targets."""
    by_tf: dict[str, set[str]] = {}
    for r, t in regulon.edges:
        if r in measured and t in measured:
            by_tf.setdefault(r, set()).add(t)

    out = []
    for tf, targets in sorted(by_tf.items(), key=lambda kv: -len(kv[1])):
        usable = sorted(targets)
        # a module needs perturbed genes, since a fold is a held-out perturbation
        n_pert_in = len((set(usable) | {tf}) & perturbed)
        if len(usable) < min_targets or n_pert_in < 2:
            continue
        genes = [tf] + usable
        if len(genes) > max_genes:
            # keep the perturbed targets first, they are what the folds can exercise
            ranked = sorted(usable, key=lambda g: (g not in perturbed, g))
            genes = [tf] + ranked[: max_genes - 1]
        out.append(Candidate(name=f"regulon_{tf}", source="regulon",
                             genes=sorted(set(genes)), seed_gene=tf))
        if len(out) >= max_modules:
            break
    return out


def coresponse_candidates(genes: list[str], responses: np.ndarray, *,
                          size: int = 20, n_modules: int = 20,
                          seed: int = 0) -> list[Candidate]:
    """Cluster genes by the correlation of their response profiles across perturbations.

    Seeds are chosen by how strongly a gene co-varies with the rest of the module, not
    by its own variance. Seeding on variance picks whichever gene is noisiest, which is
    the opposite of what a co-response module is: on a planted test where six genes
    share one program, variance-seeding built its cluster around an independent gene
    because that gene's standard deviation happened to be the largest. Summed absolute
    correlation to the other genes measures participation in a shared program directly.

    Deliberately simple otherwise: the point is a mechanical, reproducible rule, not a
    clustering contribution.
    """
    x = np.asarray(responses, float)
    if x.shape[0] < 3 or x.shape[1] < size:
        return []
    centred = x - x.mean(axis=0, keepdims=True)
    # ddof matched to the divisor below, so the diagonal is exactly 1 and the entries
    # are correlations rather than correlations scaled by n/(n-1)
    sd = centred.std(axis=0, ddof=1)
    ok = sd > 0
    corr = np.zeros((x.shape[1], x.shape[1]))
    if ok.sum() >= 2:
        sub = centred[:, ok] / sd[ok]
        c = (sub.T @ sub) / max(1, sub.shape[0] - 1)
        idx = np.flatnonzero(ok)
        corr[np.ix_(idx, idx)] = c

    connectivity = np.abs(corr).sum(axis=1) - np.abs(np.diag(corr))
    order = np.argsort(-connectivity)
    used: set[int] = set()
    out = []
    for j in order:
        if len(out) >= n_modules:
            break
        if j in used or sd[j] <= 0:
            continue
        partners = np.argsort(-np.abs(corr[j]))
        members = [int(j)]
        for p in partners:
            if len(members) >= size:
                break
            if int(p) != int(j) and int(p) not in used:
                members.append(int(p))
        used.update(members)
        out.append(Candidate(name=f"coresponse_{genes[j]}", source="coresponse",
                             genes=sorted({genes[m] for m in members}),
                             seed_gene=genes[j]))
    return out


def summarise(candidates: list[Candidate]) -> dict:
    """The reporting section 8 requires: everything extracted, pass or fail."""
    passed = [c for c in candidates if c.passed]
    reasons: dict[str, int] = {}
    for c in candidates:
        for r in c.reasons:
            key = r.split(",")[0].split(" need")[0]
            key = key.split(" ", 2)[-1] if key.startswith("only") else key
            reasons[key] = reasons.get(key, 0) + 1
    return {
        "n_extracted": len(candidates),
        "n_passed": len(passed),
        "n_failed": len(candidates) - len(passed),
        "by_source": {s: sum(1 for c in candidates if c.source == s)
                      for s in {c.source for c in candidates}},
        "passed_by_source": {s: sum(1 for c in passed if c.source == s)
                             for s in {c.source for c in candidates}},
        "failure_reasons": reasons,
    }
