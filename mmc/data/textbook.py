"""Canonical structures for the Step 1 S2 arm (PREREG_v4 section 2.2).

Hand-specified from the standard immunology of each module, written down before the
comparator runs so that S2 is a statement of textbook belief rather than a structure
tuned until it scored well. These are the edges a domain expert would draw on a
whiteboard: master transcription factors onto their signature cytokines, the
GATA3 and TBX21 mutual antagonism, the proximal TCR cascade in its accepted order,
and the CD28 and CBM route into interleukin-2.

Two limits are deliberate. The grammar allows at most three regulators per target,
so where textbook knowledge names more than three the three best-established are
kept and the rest dropped, which is recorded per module below. And an edge here
means only that the regulator's level changes the target's level, which for the
signalling genes is a statement about the transcriptional consequence of removing
the protein rather than about direct transcriptional control.

No edge in this file was added, removed, or re-signed after seeing a comparator
result. Changes require a dated amendment to the pre-registration.
"""
from __future__ import annotations

from ..grammar.model_spec import Edge, ModelSpec, Rule, Term

# (regulator, target, sign). Sign is +1 where losing the regulator lowers the target.
TEXTBOOK_EDGES: dict[str, list[tuple[str, str, int]]] = {
    # Th1 / Th2 / Th17 / Treg lineage circuit. The mutual GATA3 and TBX21 repression
    # is the load-bearing piece; the rest is master-regulator to signature cytokine.
    # IFNG in textbook terms also takes STAT4 and BATF, dropped for the in-degree cap.
    "CD4_lineage_TFs": [
        ("STAT6", "GATA3", 1),
        ("TBX21", "GATA3", -1),
        ("STAT4", "TBX21", 1),
        ("STAT1", "TBX21", 1),
        ("GATA3", "TBX21", -1),
        ("STAT3", "RORC", 1),
        ("IKZF2", "FOXP3", 1),
        ("FOXO1", "FOXP3", 1),
        ("BCL6", "PRDM1", -1),
        ("FOXO1", "TCF7", 1),
        ("TCF7", "LEF1", 1),
        ("GATA3", "IL4", 1),
        ("STAT6", "IL4", 1),
        ("IRF4", "IL4", 1),
        ("GATA3", "IL5", 1),
        ("GATA3", "IL13", 1),
        ("STAT6", "IL13", 1),
        ("TBX21", "IFNG", 1),
        ("RUNX3", "IFNG", 1),
        ("FOXP3", "IFNG", -1),
        ("RORC", "IL17A", 1),
        ("BATF", "IL17A", 1),
        ("STAT3", "IL17A", 1),
        ("RORC", "IL17F", 1),
        ("AHR", "IL17F", 1),
        ("FOXP3", "IL2", -1),
        ("PRDM1", "IL2", -1),
        ("MAF", "IL10", 1),
        ("ID2", "IL10", -1),
        ("BATF", "IL21", 1),
        ("MAF", "IL21", 1),
        ("BCL6", "IL21", 1),
        ("IRF4", "IL9", 1),
        ("GATA3", "IL9", 1),
        ("TBX21", "TNF", 1),
        ("BATF", "CSF2", 1),
        ("RORC", "CSF2", 1),
        ("FOXP3", "CTLA4", 1),
    ],
    # The proximal TCR cascade in its accepted order, terminating on interleukin-2
    # through NF-kB and AP-1. This is the module where textbook structure is least
    # ambiguous, and also the one whose leading-PC fraction is highest.
    "TCR_signalosome": [
        ("CD3E", "ZAP70", 1),
        ("ZAP70", "LAT", 1),
        ("LAT", "LCP2", 1),
        ("LAT", "PLCG1", 1),
        ("LCP2", "PLCG1", 1),
        ("PLCG1", "PRKCQ", 1),
        ("PRKCQ", "NFKB1", 1),
        ("PRKCQ", "RELA", 1),
        ("NFKB1", "RELA", 1),
        ("PRKCQ", "FOS", 1),
        ("PRKCQ", "JUN", 1),
        ("NFKB1", "IL2", 1),
        ("RELA", "IL2", 1),
        ("JUN", "IL2", 1),
    ],
    # Costimulation and the CBM complex onto the effector cytokines, plus the T-bet
    # arm. The dark candidates carry no textbook edges by construction, which is the
    # point of the arm: S2 is what was known before this atlas.
    "Cytokine_production": [
        ("CARMIL2", "CD28", 1),
        ("CD28", "IL2", 1),
        ("MALT1", "IL2", 1),
        ("IL2RB", "IL2", -1),
        ("TBX21", "IFNG", 1),
        ("MALT1", "IFNG", 1),
        ("IL2", "IFNG", 1),
        ("CD28", "CSF2", 1),
        ("MALT1", "CSF2", 1),
        ("BATF", "CSF2", 1),
        ("MALT1", "CCL3", 1),
        ("TBX21", "CCL3", 1),
        ("MALT1", "CCL4", 1),
        ("TBX21", "CCL4", 1),
        ("TBX21", "IL5", -1),
        ("BATF", "IL22", 1),
    ],
    "Th2_GATA3": [
        ("STAT6", "GATA3", 1),
        ("TBX21", "GATA3", -1),
        ("STAT4", "TBX21", 1),
        ("GATA3", "TBX21", -1),
        ("GATA3", "IL4", 1),
        ("STAT6", "IL4", 1),
        ("GATA3", "IL5", 1),
        ("GATA3", "IL13", 1),
        ("STAT6", "IL13", 1),
    ],
}


def textbook_spec(module: str, genes: list[str]) -> ModelSpec:
    """The S2 structure for a module, restricted to the genes actually in it.

    Edges naming a gene the module does not carry are dropped rather than silently
    renamed, and the count of dropped edges is available through `coverage` so the
    arm can report how much of the textbook the module could hold.
    """
    if module not in TEXTBOOK_EDGES:
        raise KeyError(f"no textbook structure recorded for {module}")
    present = set(genes)
    by_target: dict[str, list[str]] = {}
    edges: list[Edge] = []
    for reg, tgt, sign in TEXTBOOK_EDGES[module]:
        if reg not in present or tgt not in present or reg == tgt:
            continue
        edges.append(Edge(regulator=reg, target=tgt, sign=sign))
        by_target.setdefault(tgt, []).append(reg)
    rules = {t: Rule(terms=[Term(regulators=rs)]) for t, rs in by_target.items()}
    return ModelSpec(genes=list(genes), edges=edges, rules=rules)


def coverage(module: str, genes: list[str]) -> dict:
    """How much of the recorded textbook structure the module's gene set can hold."""
    total = TEXTBOOK_EDGES.get(module, [])
    present = set(genes)
    kept = [e for e in total if e[0] in present and e[1] in present and e[0] != e[1]]
    return {"module": module, "n_textbook_edges": len(total), "n_edges_in_module": len(kept),
            "dropped": [f"{r}->{t}" for r, t, _ in total if (r, t) not in
                        {(a, b) for a, b, _ in kept}]}
