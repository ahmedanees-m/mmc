"""Fetch the published Norman 2019 genetic-interaction subtype labels.

PREREG_v4 amendment A15. Step 4's pre-registered primary hypothesis names epistasis and
suppression, and amendment A7 substituted subtypes derived from this project's own fitted
coefficients because the published table is not in the GEO matrix. That substitution cannot
test a hypothesis stated about the published labels, so Step 4 is in progress until these
are obtained.

The GEO deposit carries guide identity and coverage only, confirmed by inspection, so the
labels come from the dataset loader that ships them. Only the per-perturbation subtype is
kept; the expression matrix is discarded.
"""
import json
import sys

OUT = "/work/norman_published_subtypes.json"


def main() -> int:
    import pertpy as pt

    print("loading norman_2019 through pertpy", flush=True)
    adata = pt.data.norman_2019()
    obs = adata.obs
    print("obs columns:", list(obs.columns), flush=True)

    # the subtype column has been named differently across releases, so it is found rather
    # than assumed, and what was found is recorded
    candidates = [c for c in obs.columns
                  if any(t in c.lower() for t in
                         ("genetic_interaction", "gi_", "subtype", "class", "category"))]
    print("candidate label columns:", candidates, flush=True)
    if not candidates:
        print("no subtype column present; nothing to extract", flush=True)
        return 2

    pert_col = next((c for c in obs.columns
                     if c.lower() in ("perturbation_name", "perturbation", "guide_identity",
                                      "condition", "gene_target")), None)
    print("perturbation column:", pert_col, flush=True)
    if pert_col is None:
        return 3

    out = {"source": "pertpy norman_2019", "label_column": None, "counts": {}, "labels": {}}
    for col in candidates:
        pairs = (obs[[pert_col, col]].dropna().astype(str)
                 .drop_duplicates().values.tolist())
        vals = {v for _, v in pairs}
        # the real subtype column carries the published category names
        if any(t in " ".join(vals).lower()
               for t in ("synerg", "suppress", "redund", "neomorph", "epistas")):
            out["label_column"] = col
            for p, v in pairs:
                out["labels"][p] = v
            for v in out["labels"].values():
                out["counts"][v] = out["counts"].get(v, 0) + 1
            break

    if out["label_column"] is None:
        print("candidate columns held no published subtype names:", candidates, flush=True)
        return 4

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print("column used:", out["label_column"], flush=True)
    print("counts:", out["counts"], flush=True)
    print("labelled perturbations:", len(out["labels"]), flush=True)
    print("wrote", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
