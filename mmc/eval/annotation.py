"""Scoring edge sets against annotated regulons (PREREG_v4 section 10).

Step 9 asks a different question from the rest of Step 1: not whether a structure
predicts, but whether its edges are the ones the literature already records. The two
axes come apart, and the pre-registration expects them to. An edge set can agree with
annotation and predict nothing, which is the pattern the existing record calls grounded
but non-compositional, and the oracle can predict best while agreeing least, which would
say that what predicts is not what is annotated.

Coverage is reported before any score. A precision of zero means something quite
different when the annotation has no entry for the module's genes than when it has many
and the edge set missed them, and section 10 requires the distinction be stated rather
than absorbed into the number.

CollecTRI is read from the OmniPath export, which carries gene symbols and separate
stimulation and inhibition flags, so signed and unsigned agreement can both be reported.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass


@dataclass(frozen=True)
class Regulon:
    """An annotated regulatory network, restricted to a gene universe on demand."""

    edges: set[tuple[str, str]]
    signed: dict[tuple[str, str], int]
    name: str = "collectri"

    def restrict(self, genes) -> Regulon:
        universe = set(genes)
        keep = {(r, t) for r, t in self.edges if r in universe and t in universe}
        return Regulon(keep, {k: v for k, v in self.signed.items() if k in keep},
                       self.name)

    def coverage(self, genes) -> dict:
        """How much of the module the annotation actually speaks to."""
        universe = set(genes)
        regs = {r for r, _ in self.edges} & universe
        tgts = {t for _, t in self.edges} & universe
        within = self.restrict(universe)
        return {
            "n_module_genes": len(universe),
            "n_genes_as_annotated_regulator": len(regs),
            "n_genes_as_annotated_target": len(tgts),
            "n_annotated_edges_within_module": len(within.edges),
            "n_signed_within_module": len(within.signed),
            "fraction_genes_covered": (len(regs | tgts) / len(universe)) if universe else 0.0,
        }


def load_collectri(path: str) -> Regulon:
    """Read the OmniPath CollecTRI export.

    Sign is taken only where the annotation is unambiguous: an interaction flagged as
    both stimulating and inhibiting, or as neither, carries no sign and is counted in
    the unsigned edge set alone.
    """
    edges: set[tuple[str, str]] = set()
    signed: dict[tuple[str, str], int] = {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            r = (row.get("source_genesymbol") or "").strip()
            t = (row.get("target_genesymbol") or "").strip()
            if not r or not t or r == t:
                continue
            edges.add((r, t))
            stim = (row.get("is_stimulation") or "").strip().lower() in ("true", "1")
            inhib = (row.get("is_inhibition") or "").strip().lower() in ("true", "1")
            if stim and not inhib:
                signed[(r, t)] = 1
            elif inhib and not stim:
                signed[(r, t)] = -1
    return Regulon(edges, signed)


def score_edges(edge_set, regulon: Regulon, genes, *, use_sign: bool = False) -> dict:
    """Precision, recall and Jaccard of an edge set against the annotation.

    The annotation is restricted to the module's genes first, so recall is against what
    could have been recovered rather than against the whole regulon, which no module
    could ever match.
    """
    within = regulon.restrict(genes)
    truth = set(within.signed) if use_sign else set(within.edges)
    if use_sign:
        proposed = {(r, t) for r, t, s in edge_set if within.signed.get((r, t)) == s}
        proposed_all = {(r, t) for r, t, _ in edge_set}
        hit = proposed & truth
        n_proposed = len(proposed_all)
    else:
        proposed_all = {(r, t) for r, t, _ in edge_set}
        hit = proposed_all & truth
        n_proposed = len(proposed_all)

    precision = len(hit) / n_proposed if n_proposed else float("nan")
    recall = len(hit) / len(truth) if truth else float("nan")
    union = len(proposed_all | truth)
    return {
        "n_proposed": n_proposed,
        "n_annotated_within_module": len(truth),
        "n_hit": len(hit),
        "precision": precision,
        "recall": recall,
        "jaccard": (len(hit) / union) if union else float("nan"),
        "signed": use_sign,
    }


def enrichment(edge_set, regulon: Regulon, genes, *, n_draws: int = 2000,
               seed: int = 0) -> dict:
    """Is the agreement better than a random edge set of the same size would give?

    Precision against a sparse annotation is hard to read on its own: a module whose
    annotation holds three edges out of a possible 756 makes a precision of 0.05 look
    poor when it is many times chance. This reports the chance level explicitly.
    """
    import numpy as np

    within = regulon.restrict(genes)
    truth = set(within.edges)
    universe = [(r, t) for r in genes for t in genes if r != t]
    if not truth or not universe or not edge_set:
        return {"chance_precision": float("nan"), "ratio": float("nan"), "p": float("nan")}

    observed = len({(r, t) for r, t, _ in edge_set} & truth)
    k = len({(r, t) for r, t, _ in edge_set})
    rng = np.random.default_rng(seed)
    idx = np.arange(len(universe))
    hits = []
    for _ in range(n_draws):
        pick = rng.choice(idx, size=min(k, len(universe)), replace=False)
        hits.append(len({universe[i] for i in pick} & truth))
    hits = np.asarray(hits, float)
    chance = float(hits.mean()) / k if k else float("nan")
    obs_prec = observed / k if k else float("nan")
    return {
        "observed_hits": observed,
        "chance_hits": float(hits.mean()),
        "chance_precision": chance,
        "observed_precision": obs_prec,
        "ratio": (obs_prec / chance) if chance else float("inf") if observed else float("nan"),
        "p": float(((hits >= observed).sum() + 1) / (n_draws + 1)),
    }
