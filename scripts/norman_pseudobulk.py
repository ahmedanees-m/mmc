"""Work Stream D, step 1: pseudobulk the Norman 2019 CRISPRa Perturb-seq matrix.

Streams the GSE133344 filtered 10x matrix (33694 genes x 111668 cells, ~362M nonzeros)
in a single low-memory pass, groups cells by perturbation identity, and writes a
perturbation-by-gene log-fold-change table (versus the pooled non-targeting control) for
the downstream epistasis compose test.

Perturbation identity (cell_identities guide_identity, form "A_B__A_B"):
  both tokens NegCtrl        -> control (pooled reference)
  one gene + one NegCtrl     -> single perturbation of that gene (both guide orders pooled)
  two genes                  -> double perturbation of that unordered pair

Only cells with good_coverage and a single confident identity are used. Reads from /norman,
writes /norman/norman_pseudobulk.npz.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

NORM = "/norman"
GENES = f"{NORM}/GSE133344_filtered_genes.tsv.gz"
BARCODES = f"{NORM}/GSE133344_filtered_barcodes.tsv.gz"
IDENT = f"{NORM}/GSE133344_filtered_cell_identities.csv.gz"
MTX = f"{NORM}/GSE133344_filtered_matrix.mtx.gz"
OUT = f"{NORM}/norman_pseudobulk.npz"

LFC_MIN = 0.5          # a gene is a readout if |logFC| >= this in any non-control group
CHUNK = 20_000_000     # matrix nonzeros per streamed chunk


def classify(gid: str):
    """Return (type, group_label, geneA, geneB) or None for an unparseable label."""
    left = gid.split("__", 1)[0]
    toks = left.split("_")
    if len(toks) != 2:
        return None
    a, b = toks
    a_ctrl, b_ctrl = a.startswith("NegCtrl"), b.startswith("NegCtrl")
    if a_ctrl and b_ctrl:
        return ("ctrl", "CTRL", "", "")
    if a_ctrl:
        return ("single", f"S:{b}", b, "")
    if b_ctrl:
        return ("single", f"S:{a}", a, "")
    x, y = sorted([a, b])
    if x == y:
        return None
    return ("double", f"D:{x}+{y}", x, y)


def main() -> None:
    gene_names = pd.read_csv(GENES, sep="\t", header=None)[1].astype(str).values
    n_genes = len(gene_names)
    barcodes = pd.read_csv(BARCODES, sep="\t", header=None)[0].astype(str).values
    n_cells = len(barcodes)
    bc_to_col = {bc: i for i, bc in enumerate(barcodes)}
    print(f"genes {n_genes}  cells {n_cells}", flush=True)

    ident = pd.read_csv(IDENT)
    ident = ident[(ident["good_coverage"]) & (ident["number_of_cells"] == 1)].copy()
    parsed = ident["guide_identity"].map(classify)
    bad = int(parsed.isna().sum())
    ident = ident[parsed.notna()].copy()
    parsed = parsed[parsed.notna()]
    ident["gtype"] = [p[0] for p in parsed]
    ident["glabel"] = [p[1] for p in parsed]
    ident["gA"] = [p[2] for p in parsed]
    ident["gB"] = [p[3] for p in parsed]
    print(f"clean cells {len(ident)}  unparseable labels {bad}", flush=True)

    groups = sorted(ident["glabel"].unique().tolist())
    gidx = {g: i for i, g in enumerate(groups)}
    n_grp = len(groups)
    meta = ident.drop_duplicates("glabel").set_index("glabel")
    gtype = np.array([meta.loc[g, "gtype"] for g in groups], dtype="U8")
    gA = np.array([meta.loc[g, "gA"] for g in groups], dtype="U20")
    gB = np.array([meta.loc[g, "gB"] for g in groups], dtype="U20")
    print(f"groups {n_grp}  (singles {int((gtype=='single').sum())}, "
          f"doubles {int((gtype=='double').sum())}, ctrl {int((gtype=='ctrl').sum())})",
          flush=True)

    col_group = np.full(n_cells, -1, dtype=np.int64)
    ncells = np.zeros(n_grp, dtype=np.int64)
    for bc, gl in zip(ident["cell_barcode"].values, ident["glabel"].values):
        c = bc_to_col.get(bc)
        if c is not None:
            col_group[c] = gidx[gl]
            ncells[gidx[gl]] += 1

    acc = np.zeros(n_genes * n_grp, dtype=np.float64)
    reader = pd.read_csv(MTX, sep=" ", skiprows=3, header=None, names=["g", "c", "v"],
                         dtype={"g": np.int64, "c": np.int64, "v": np.float64},
                         chunksize=CHUNK)
    total = 0
    for k, chunk in enumerate(reader):
        gi = chunk["g"].values - 1
        ci = chunk["c"].values - 1
        v = chunk["v"].values
        if k == 0:
            assert not np.isnan(v).any(), "matrix parse produced NaN (delimiter?)"
            assert gi.max() < n_genes and ci.max() < n_cells, "index out of range"
        grp = col_group[ci]
        keep = grp >= 0
        flat = gi[keep] * n_grp + grp[keep]
        acc += np.bincount(flat, weights=v[keep], minlength=n_genes * n_grp)
        total += len(chunk)
        print(f"  chunk {k}  rows {total}", flush=True)
    acc = acc.reshape(n_genes, n_grp)
    print(f"streamed {total} nonzeros", flush=True)

    tot = acc.sum(axis=0)
    cpm = np.divide(acc, tot[None, :], out=np.zeros_like(acc), where=tot[None, :] > 0) * 1e6
    logexpr = np.log2(1.0 + cpm)
    ctrl = logexpr[:, gidx["CTRL"]]
    logfc = logexpr - ctrl[:, None]                      # genes x groups

    nonctrl = gtype != "ctrl"
    readout_mask = (np.abs(logfc[:, nonctrl]).max(axis=1) >= LFC_MIN)
    ro = np.where(readout_mask)[0]
    print(f"readout genes (|logFC|>={LFC_MIN} in any group) {len(ro)}", flush=True)

    np.savez_compressed(
        OUT,
        readout_genes=gene_names[ro].astype("U20"),
        group_label=np.array(groups, dtype="U40"),
        group_type=gtype, gene_A=gA, gene_B=gB,
        n_cells=ncells,
        logfc=logfc[ro].T.astype(np.float32),            # groups x readout
    )
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
