"""Accessor over the existing VERDICT Zhu differential-expression store.

MMC does not re-ingest the atlas. It queries the DuckDB and Parquet store already
built for VERDICT, pointed at by MMC_ZHU_STORE (the directory holding
perturbation.parquet, gene.parquet, and de.parquet).

Store schema:
    perturbation.parquet  pert_code, perturbation, perturbation_id, condition,
                          n_cells_target, guide_correlation_signif, n_downstream,
                          n_downstream_pct, ...
    gene.parquet          gene_code, gene_symbol, gene_ensembl
    de.parquet            pert_code, gene_code, log_fc, adj_p_value, zscore

A knockdown effect is de.log_fc for the pair (pert_code = regulator knockdown,
gene_code = target); adj_p_value is the FDR; guide_correlation_signif is the
perturbation's cross-guide concordance (often null for single-guide perturbations).
Symbols are resolved to integer codes first so the large de table is filtered on
its join keys rather than scanned.
"""
from __future__ import annotations

import os
from functools import lru_cache

import duckdb
import pandas as pd

CONDITIONS = ("Rest", "Stim8hr", "Stim48hr")
_EFFECT_COLS = ["perturbation", "target_gene", "effect_size", "fdr", "crossguide_r"]


def store_path() -> str:
    p = os.environ.get("MMC_ZHU_STORE")
    if not p:
        raise RuntimeError(
            "Set MMC_ZHU_STORE to the existing Zhu store directory "
            "(holding perturbation.parquet, gene.parquet, de.parquet)."
        )
    return p


@lru_cache(maxsize=1)
def _con() -> duckdb.DuckDBPyConnection:
    base = store_path().rstrip("/")
    con = duckdb.connect()
    con.execute(f"CREATE VIEW zhu_pert AS SELECT * FROM read_parquet('{base}/perturbation.parquet')")
    con.execute(f"CREATE VIEW zhu_gene AS SELECT * FROM read_parquet('{base}/gene.parquet')")
    con.execute(f"CREATE VIEW zhu_de   AS SELECT * FROM read_parquet('{base}/de.parquet')")
    return con


def _pert_codes(regulators: list[str], condition: str) -> dict[int, str]:
    """Map pert_code to regulator symbol for the perturbations in this condition."""
    regs = [g for g in regulators if g]
    if not regs:
        return {}
    ph = ",".join("?" * len(regs))
    rows = _con().execute(
        f"SELECT perturbation, pert_code FROM zhu_pert "
        f"WHERE condition = ? AND (perturbation IN ({ph}) OR perturbation_id IN ({ph}))",
        [condition] + regs + regs,
    ).fetchall()
    return {int(c): str(s) for s, c in rows}


def _gene_codes(targets: list[str]) -> dict[int, str]:
    """Map gene_code to target symbol for the requested genes."""
    tgts = [g for g in targets if g]
    if not tgts:
        return {}
    ph = ",".join("?" * len(tgts))
    rows = _con().execute(
        f"SELECT gene_code, gene_symbol FROM zhu_gene "
        f"WHERE gene_symbol IN ({ph}) OR gene_ensembl IN ({ph})",
        tgts + tgts,
    ).fetchall()
    return {int(c): str(s) for c, s in rows}


def module_effects(regulators: list[str], targets: list[str],
                   condition: str) -> pd.DataFrame:
    """Signed knockdown effect of each regulator on each target, one condition.

    Columns: perturbation, target_gene, effect_size (log2 fold change), fdr,
    crossguide_r (the perturbation's cross-guide concordance). This is the atom of
    both the precondition test and the Step-1 module slice.
    """
    assert condition in CONDITIONS, condition
    pcodes = _pert_codes(regulators, condition)
    gcodes = _gene_codes(targets)
    if not pcodes or not gcodes:
        return pd.DataFrame(columns=_EFFECT_COLS)

    pph = ",".join(str(c) for c in pcodes)
    gph = ",".join(str(c) for c in gcodes)
    de = _con().execute(
        f"SELECT pert_code, gene_code, log_fc, adj_p_value FROM zhu_de "
        f"WHERE pert_code IN ({pph}) AND gene_code IN ({gph})"
    ).df()
    xg = _con().execute(
        f"SELECT pert_code, guide_correlation_signif FROM zhu_pert "
        f"WHERE pert_code IN ({pph})"
    ).df()
    xg_map = dict(zip(xg["pert_code"].astype(int), xg["guide_correlation_signif"]))

    de["perturbation"] = de["pert_code"].astype(int).map(pcodes)
    de["target_gene"] = de["gene_code"].astype(int).map(gcodes)
    de["crossguide_r"] = de["pert_code"].astype(int).map(xg_map)
    de = de.rename(columns={"log_fc": "effect_size", "adj_p_value": "fdr"})
    return de[_EFFECT_COLS]


def perturbation_response(gene: str, condition: str) -> pd.DataFrame:
    """All measured-gene deltas for one knockdown in one condition, for fitting or eval.

    Columns: target_gene, effect_size (log2 fold change), fdr, crossguide_r, n_cells.
    """
    assert condition in CONDITIONS, condition
    row = _con().execute(
        "SELECT pert_code, guide_correlation_signif, n_cells_target FROM zhu_pert "
        "WHERE condition = ? AND (perturbation = ? OR perturbation_id = ?) LIMIT 1",
        [condition, gene, gene],
    ).fetchone()
    if row is None:
        return pd.DataFrame(columns=["target_gene", "effect_size", "fdr", "crossguide_r", "n_cells"])
    pert_code, crossguide_r, n_cells = int(row[0]), row[1], row[2]
    df = _con().execute(
        "SELECT g.gene_symbol AS target_gene, d.log_fc AS effect_size, "
        "       d.adj_p_value AS fdr "
        "FROM zhu_de d JOIN zhu_gene g ON d.gene_code = g.gene_code "
        "WHERE d.pert_code = ?",
        [pert_code],
    ).df()
    df["crossguide_r"] = crossguide_r
    df["n_cells"] = n_cells
    return df


def measured_genes() -> set[str]:
    """The set of gene symbols measured in the atlas, for coverage checks."""
    rows = _con().execute("SELECT gene_symbol FROM zhu_gene").fetchall()
    return {str(r[0]) for r in rows}
