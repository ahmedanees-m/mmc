"""Amendment A19: the identifiability signature on data this project did not generate.

Diagnostics only. No model is fitted, no comparator is run, and no hypothesis is tested
about these datasets. The question is narrow: is the low effective rank that explains the
result on the primary module a property of that module set, or a property of perturbation
response matrices more generally.

Two sources.

  zhu      the remaining states of the atlas used throughout, Rest and Stim48hr, on the
           same modules, so the only thing that changes is the cell state
  norman   single perturbations from the Norman combinatorial CRISPRa atlas, already in
           pseudobulk from the work in section 6. A different cell type, a different
           laboratory, and gain of function rather than loss of function

A dataset that does not show the signature is as informative as one that does, and is
reported without adjustment.

    python scripts/step13_external_diagnostics.py --source zhu \\
        --module coresponse_PIM1 --module-def /work/step7_defs/coresponse_PIM1.json \\
        --out /work/results/a19/zhu_coresponse_PIM1.json
    python scripts/step13_external_diagnostics.py --source norman \\
        --npz /n/norman_pseudobulk.npz --out /work/results/a19/norman_singles.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mmc.eval import identifiability as ident


def diagnostics_of(matrix: np.ndarray) -> dict:
    m = np.asarray(matrix, float)
    return {
        "n_perturbations": int(m.shape[0]),
        "n_genes": int(m.shape[1]),
        "effective_rank": float(ident.effective_rank(m)),
        "effective_rank_normalised": float(ident.effective_rank(m) / max(1, m.shape[0])),
        "leading_pc_fraction": float(ident.mean_leading_pc_fraction(m)),
        "perturbation_specific_ratio": float(ident.perturbation_specific_ratio(m)),
    }


def from_zhu(args) -> dict:
    from mmc.data import module_data, module_extract

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(d["module"], d["regulators"], d["targets"])

    out = {"source": "zhu", "module": args.module, "states": {}}
    for condition in args.conditions.split(","):
        try:
            mod = module_data.build_module_data(args.module, condition)
        except Exception as e:                       # a module may not be measured in a state
            out["states"][condition] = {"error": str(e)}
            print(f"  {condition}: unavailable, {e}")
            continue
        d = diagnostics_of(np.asarray(mod.observed, float))
        out["states"][condition] = d
        print(f"  {condition}: {d['n_perturbations']} perturbations, "
              f"effective rank {d['effective_rank']:.3f} "
              f"({d['effective_rank_normalised']:.3f} normalised), "
              f"leading PC {d['leading_pc_fraction']:.4f}")
    return out


def matched_submatrices(z, args) -> dict:
    """Norman diagnostics on submatrices shaped like the modules used here.

    Added 2026-08-22, after the whole-transcriptome numbers came back. The modules in this
    study are 11 to 40 gene sets whose readouts are the perturbed genes themselves, while
    the Norman matrix is 105 perturbations by 20,421 genes. Effective rank is not
    comparable across those shapes, so a difference between them says nothing on its own.
    This draws perturbation-by-gene submatrices of the same shape as the modules, with the
    readout set equal to the perturbed set exactly as `module_extract` builds them, and
    reports the distribution over draws.
    """
    gt = np.asarray(z["group_type"])
    logfc = np.asarray(z["logfc"], float)
    readout = np.asarray(z["readout_genes"]).astype(str)
    gene_a = np.asarray(z["gene_A"]).astype(str)

    singles = gt == "single"
    perturbed = gene_a[singles]
    rows = logfc[singles]
    index = {g: i for i, g in enumerate(readout)}
    usable = [(i, g) for i, g in enumerate(perturbed) if g in index]

    rng = np.random.default_rng(0)
    out = {}
    for k in (11, 20, 28, 40):
        if k > len(usable):
            continue
        draws = []
        for _ in range(args.draws):
            pick = rng.choice(len(usable), size=k, replace=False)
            row_idx = [usable[i][0] for i in pick]
            col_idx = [index[usable[i][1]] for i in pick]
            draws.append(diagnostics_of(rows[np.ix_(row_idx, col_idx)]))
        out[f"k={k}"] = {
            "n_draws": len(draws),
            "effective_rank_mean": float(np.mean([d["effective_rank"] for d in draws])),
            "effective_rank_sd": float(np.std([d["effective_rank"] for d in draws])),
            "effective_rank_normalised_mean":
                float(np.mean([d["effective_rank_normalised"] for d in draws])),
            "leading_pc_fraction_mean":
                float(np.nanmean([d["leading_pc_fraction"] for d in draws])),
        }
        r = out[f"k={k}"]
        print(f"    {k} by {k}: effective rank {r['effective_rank_mean']:.3f} "
              f"(sd {r['effective_rank_sd']:.3f}, "
              f"{r['effective_rank_normalised_mean']:.3f} normalised), "
              f"leading PC {r['leading_pc_fraction_mean']:.4f}")
    return out


def from_norman(args) -> dict:
    z = np.load(args.npz, allow_pickle=True)
    gt = np.asarray(z["group_type"])
    logfc = np.asarray(z["logfc"], float)

    out = {"source": "norman", "npz": args.npz, "subsets": {}}
    for label, mask in (("single", gt == "single"), ("double", gt == "double")):
        m = logfc[mask]
        if m.shape[0] < 3:
            continue
        # the full transcriptome, then restricted to the genes that actually move, since a
        # matrix dominated by unexpressed genes would flatter the diagnostic
        d_all = diagnostics_of(m)
        moving = np.abs(m).max(axis=0) > 0.5
        d_de = diagnostics_of(m[:, moving])
        out["subsets"][label] = {"all_genes": d_all, "responsive_genes": d_de,
                                 "n_responsive": int(moving.sum())}
        print(f"  {label}: {m.shape[0]} perturbations")
        print(f"    all {d_all['n_genes']} genes:        effective rank "
              f"{d_all['effective_rank']:.3f} ({d_all['effective_rank_normalised']:.3f} "
              f"normalised), leading PC {d_all['leading_pc_fraction']:.4f}")
        print(f"    {int(moving.sum())} responsive genes: effective rank "
              f"{d_de['effective_rank']:.3f} ({d_de['effective_rank_normalised']:.3f} "
              f"normalised), leading PC {d_de['leading_pc_fraction']:.4f}")

    print("  shape-matched submatrices, readouts equal to the perturbed set:")
    out["shape_matched"] = matched_submatrices(z, args)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=("zhu", "norman"), required=True)
    ap.add_argument("--module", default="")
    ap.add_argument("--module-def", default="")
    ap.add_argument("--conditions", default="Rest,Stim8hr,Stim48hr")
    ap.add_argument("--npz", default="")
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print(f"=== A19 diagnostics, source {args.source} "
          f"{args.module or Path(args.npz).name} ===", flush=True)
    out = from_zhu(args) if args.source == "zhu" else from_norman(args)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
