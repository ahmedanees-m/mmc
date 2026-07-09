"""Spec the powered cytokine-production module from the Zhu store, and split its
regulators into the KNOWN backbone vs the DARK discovery candidates.

Per v3 Part 3a / Part 4: assemble a module that (a) has enough DE power for the held-out
gate and edge-ablation to discriminate, and (b) deliberately includes data-supported
regulators that are NOT canonical and NOT named by Zhu, because a discovery can only come
from that dark set. A module of only known regulators has power and a clean fit but
cannot yield a novel discovery.

ZHU_NAMED_REGULATORS is completed from the Zhu 2025 paper and its analysis repository
(github.com/emdann/GWT_perturbseq_analysis_2025): the abstract names the SAGA and Mediator
complexes as novel cytokine regulators, and the repository's arrayed-FACS and bulk-RNAseq
validation tables (IL10_IL21_arrayed_validation.csv, IL10IL21bulkRNAseq_DESeq2_results.csv)
name ATP2A2, CYB5R4, ELOB, MEN1, KDM1A, SGF29, MED24 as validated IL10/IL21 regulators.
The full SAGA and Mediator complexes are excluded because Zhu reported the complexes.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

# ----------------------------- panels & lists -----------------------------
# Broad cytokine/chemokine readout (the DE columns) -- breadth = power.
CYTOKINE_PANEL = [
    "IL2", "IL3", "IL4", "IL5", "IL9", "IL10", "IL13", "IL21", "IL17A", "IL17F", "IL22",
    "IFNG", "TNF", "LTA", "CSF2", "IL16", "CCL3", "CCL4", "CCL5", "CXCL8", "TGFB1",
]

# Canonical / textbook cytokine regulators -- KNOWN backbone (subtract for discovery).
CANONICAL_REGULATORS = {
    "GATA3", "TBX21", "RORC", "RORA", "FOXP3", "BCL6", "PRDM1",
    "STAT1", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "JUN", "JUNB", "FOS", "FOSB", "BATF", "BATF3", "MAF", "MAFB",
    "IRF1", "IRF4", "IRF8",
    "NFATC1", "NFATC2", "NFKB1", "NFKB2", "RELA", "RELB", "REL",
    "AHR", "RUNX1", "RUNX3", "ETS1", "EGR2", "NR4A1", "NR4A2", "NR4A3", "TOX",
    "ID2", "ID3", "MYB", "GFI1", "IKZF1", "IKZF2", "IKZF3", "POU2F2", "NFE2L2",
    "JAK1", "JAK2", "JAK3", "TYK2",
}

# Regulators NAMED by Zhu 2025 as cytokine regulators -- excluded from discovery, or a
# discovery re-derives the atlas paper. Verified from the abstract (SAGA + Mediator drive
# cytokines including TNF and IL16) and the Zhu analysis repository's validation tables
# (ATP2A2, CYB5R4, ELOB, MEN1, KDM1A, SGF29, MED24 validated for IL10/IL21).
_SAGA = {
    "KAT2A", "KAT2B", "TADA1", "TADA2A", "TADA2B", "TADA3", "TAF5L", "TAF6L", "TAF9",
    "TAF10", "TAF12", "TAF9B", "SUPT3H", "SUPT7L", "SUPT20H", "CCDC101", "SGF29",
    "ATXN7", "ATXN7L1", "ATXN7L2", "ATXN7L3", "ENY2", "USP22",
}
_MEDIATOR = {
    "MED1", "MED4", "MED6", "MED7", "MED8", "MED9", "MED10", "MED11", "MED12", "MED12L",
    "MED13", "MED13L", "MED14", "MED15", "MED16", "MED17", "MED18", "MED19", "MED20",
    "MED21", "MED22", "MED23", "MED24", "MED25", "MED26", "MED27", "MED28", "MED29",
    "MED30", "MED31", "CCNC", "CDK8", "CDK19",
}
_VALIDATED = {"ATP2A2", "CYB5R4", "ELOB", "MEN1", "KDM1A"}
ZHU_NAMED_REGULATORS = _SAGA | _MEDIATOR | _VALIDATED

# Textbook TCR signalosome, costimulation, and cytokine/interleukin receptors. Knocking
# these down lowers cytokine production because they are upstream of activation, so they
# have strong cytokine effects but are not novel; they are part of the KNOWN backbone
# (scaffold and power), never the dark discovery set. Without this, the dark list is
# dominated by the signalosome (CD3, ZAP70, LAT, LCP2, PLCG1, VAV1), which is textbook.
KNOWN_SIGNALING = {
    "CD3D", "CD3E", "CD3G", "CD247", "ZAP70", "LCK", "FYN", "LAT", "LCP2", "PLCG1",
    "PLCG2", "VAV1", "VAV2", "VAV3", "ITK", "PRKCQ", "GRAP2", "THEMIS", "RASGRP1",
    "SOS1", "SOS2", "WAS", "WIPF1", "NCK1", "NCK2", "RLTPR", "CARMIL2", "DEF6",
    "CARD11", "MALT1", "BCL10", "CD28", "ICOS", "PDPK1", "PRKCB", "MAP3K7", "CHUK",
    "IKBKB", "IKBKG", "LCP1", "FYB1", "SKAP1", "INPP5D", "PTPN6", "PTPN22", "CSK",
    "IL2RA", "IL2RB", "IL2RG", "IL4R", "IL6R", "IL6ST", "IL7R", "IL10RA", "IL10RB",
    "IL12RB1", "IL12RB2", "IL18R1", "IL18RAP", "IL21R", "IL23R", "IFNGR1", "IFNGR2",
    "TNFRSF1A", "TNFRSF1B", "IL1R1", "CSF2RB", "IL9R", "IL13RA1", "IL15RA",
}

FDR_THRESH = 0.10
MIN_DE_ENTRIES = 40     # power precondition: total significant (reg x cytokine) LOO entries
MODULE_MAX_DARK = 15    # top dark candidates to include (tractability vs discovery room)


# ----------------------------- store query --------------------------------
def query_cytokine_hits(conditions=("Rest", "Stim8hr", "Stim48hr")) -> pd.DataFrame:
    """Per (perturbation, cytokine, condition): effect_size, fdr, for perturbations that
    significantly affect a cytokine. Adapted to the raw Zhu store (codes joined to names)."""
    from mmc.shared.store import _con
    con = _con()
    cyq = ",".join("?" for _ in CYTOKINE_PANEL)
    gcodes = [int(c) for (c,) in con.execute(
        f"SELECT gene_code FROM zhu_gene WHERE gene_symbol IN ({cyq})", CYTOKINE_PANEL).fetchall()]
    if not gcodes:
        return pd.DataFrame(columns=["perturbation", "cytokine", "condition", "effect_size", "fdr"])
    gph = ",".join(str(c) for c in gcodes)
    condq = ",".join("?" for _ in conditions)
    df = con.execute(
        f"SELECT p.perturbation AS perturbation, g.gene_symbol AS cytokine, "
        f"       p.condition AS condition, d.log_fc AS effect_size, d.adj_p_value AS fdr "
        f"FROM zhu_de d "
        f"JOIN zhu_pert p ON d.pert_code = p.pert_code "
        f"JOIN zhu_gene g ON d.gene_code = g.gene_code "
        f"WHERE d.gene_code IN ({gph}) AND d.adj_p_value < {FDR_THRESH} "
        f"      AND p.condition IN ({condq})",
        list(conditions),
    ).df()
    df = df[df["perturbation"] != df["cytokine"]]                 # drop self-effects
    df = df[~df["perturbation"].isin(CYTOKINE_PANEL)]             # a readout is not a regulator
    return df


# ----------------------- partition + power + assembly ---------------------
def build_module(hits: pd.DataFrame):
    """From significant cytokine hits, rank regulators and split known vs dark."""
    sig = hits.copy()
    sig["abs_eff"] = sig["effect_size"].abs()
    g = (sig.groupby("perturbation")
            .agg(n_cyto=("cytokine", "nunique"),
                 n_entries=("cytokine", "size"),
                 mean_abs_eff=("abs_eff", "mean"))
            .reset_index()
            .sort_values(["n_cyto", "mean_abs_eff"], ascending=False))
    # Zhu-named is checked first (exclusion takes priority), then canonical TFs, then
    # textbook signalling/receptors; anything else with a cytokine effect is dark.
    conds = [g["perturbation"].isin(ZHU_NAMED_REGULATORS),
             g["perturbation"].isin(CANONICAL_REGULATORS),
             g["perturbation"].isin(KNOWN_SIGNALING)]
    g["class"] = np.select(conds, ["zhu_named", "canonical", "signaling"], default="dark")
    backbone = g[g["class"].isin(["canonical", "signaling"])]
    dark = g[g["class"] == "dark"].head(MODULE_MAX_DARK)
    module_regs = pd.concat([backbone, dark])
    total_entries = int(module_regs["n_entries"].sum())
    return g, backbone, dark, total_entries


def report(g, backbone, dark, total_entries):
    print(f"cytokine panel: {len(CYTOKINE_PANEL)} readout genes")
    print(f"regulators with significant cytokine effects: {len(g)}")
    print(f"  canonical TFs: {(g['class']=='canonical').sum()} | signalling/receptors: "
          f"{(g['class']=='signaling').sum()} | dark candidates: {(g['class']=='dark').sum()} "
          f"| zhu_named (excluded): {(g['class']=='zhu_named').sum()}")
    print(f"\nPOWER PRECONDITION: {total_entries} total DE entries (need >= {MIN_DE_ENTRIES}) "
          f"-> {'PASS' if total_entries>=MIN_DE_ENTRIES else 'FAIL: widen panel / add regulators'}")
    print("\nTop DARK discovery candidates (data-supported, not canonical, not Zhu-named):")
    print(dark[["perturbation", "n_cyto", "n_entries", "mean_abs_eff"]].to_string(index=False))
    zn = g[g["class"] == "zhu_named"]
    if len(zn):
        print("\nExcluded because Zhu named them (would be re-derivation):")
        print(zn[["perturbation", "n_cyto", "n_entries", "mean_abs_eff"]].head(15).to_string(index=False))


def main():
    hits = query_cytokine_hits()
    g, backbone, dark, total = build_module(hits)
    report(g, backbone, dark, total)
    module_genes = sorted(set(CYTOKINE_PANEL) | set(backbone["perturbation"]) | set(dark["perturbation"]))
    pd.Series(module_genes).to_csv("/app/cytokine_module_genes.csv", index=False, header=["gene"])
    print(f"\nwrote {len(module_genes)} module genes -> cytokine_module_genes.csv")


# ------------------------------- self test --------------------------------
def selftest():
    rows = []

    def add(reg, cys, eff):
        for c in cys:
            rows.append((reg, c, "Stim8hr", eff, 0.01))

    add("GATA3", ["IL5", "IL13", "IL4"], -3.0)          # canonical
    add("STAT6", ["IL4", "IL13"], -2.0)                 # canonical
    add("BATF", ["IL17A", "IL21", "IFNG"], -1.5)        # canonical
    add("TADA2B", ["TNF", "IL16"], -1.0)                # zhu_named (SAGA) -> excluded
    add("USP22", ["TNF"], -0.8)                         # zhu_named (SAGA)
    add("CYB5R4", ["IL10"], -0.9)                       # zhu_named (validated) -> excluded
    add("DARK1", ["IL5", "IL13", "IL21", "IFNG"], -2.5)  # dark, broad -> top candidate
    add("DARK2", ["IL2", "TNF"], -1.2)                  # dark
    add("DARK3", ["IL10"], -0.9)                        # dark
    add("DARK4", ["CSF2", "IL9"], -1.1)                 # dark
    hits = pd.DataFrame(rows, columns=["perturbation", "cytokine", "condition", "effect_size", "fdr"])
    g, backbone, dark, total = build_module(hits)
    report(g, backbone, dark, total)
    assert set(dark["perturbation"]).issubset({"DARK1", "DARK2", "DARK3", "DARK4"})
    assert "TADA2B" not in set(dark["perturbation"]) and "CYB5R4" not in set(dark["perturbation"])
    assert dark.iloc[0]["perturbation"] == "DARK1"
    print("\n[selftest] partition + ranking + Zhu-exclusion correct.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    selftest() if args.selftest else main()
