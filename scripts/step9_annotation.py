"""Step 9: annotation agreement against held-out predictive advantage (PREREG_v4 section 10).

Scores every structure source from Step 1 against CollecTRI and places each on two axes:
how much it agrees with the recorded literature, and how much predictive advantage it
carries. The pre-registration expects these to come apart, and which way they come apart
is the result.

Coverage is printed before any score. Precision against a sparse annotation is not
interpretable on its own, so a permutation chance level is reported alongside every
number: a precision of 0.05 is poor if chance is 0.04 and notable if chance is 0.005.

    python scripts/step9_annotation.py --results 'results/step1_*.json' \\
        --collectri /work/data/collectri.tsv --out results/step9.json
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

from mmc.eval import annotation as ann

STRUCTURE_SOURCES = ("oracle", "oracle_nested", "textbook", "claude",
                     "mean_difference", "grnboost2")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--collectri", required=True)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    regulon = ann.load_collectri(args.collectri)
    print(f"CollecTRI: {len(regulon.edges)} interactions, "
          f"{len(regulon.signed)} with an unambiguous sign")

    paths = [p for pattern in args.results for p in glob.glob(pattern)]
    out = {"collectri_edges": len(regulon.edges), "modules": {}}

    for path in sorted(paths):
        with open(path) as f:
            d = json.load(f)
        module = d.get("module")
        if not module or "specs" not in d:
            continue
        genes = d["module_summary"]["perts"]
        spec_any = next(iter(d["specs"].values()), None)
        if spec_any:
            genes = spec_any["genes"]
        cov = regulon.coverage(genes)
        entry = {"coverage": cov, "sources": {}}
        print(f"\n=== {module} ===")
        print(f"  coverage: {cov['n_annotated_edges_within_module']} annotated edges "
              f"within the module, {cov['fraction_genes_covered']:.0%} of genes appear "
              f"in CollecTRI at all")
        if cov["n_annotated_edges_within_module"] == 0:
            print("  no annotated edge lies inside this module; agreement is undefined "
                  "here and is reported as such rather than as zero")

        # predictive advantage, for the second axis
        deltas = {r["source"]: r["delta_vs_linear"] for r in d.get("table", [])}

        for name in STRUCTURE_SOURCES:
            spec = d["specs"].get(name)
            if not spec:
                continue
            edges = [(e["regulator"], e["target"], e["sign"]) for e in spec["edges"]]
            unsigned = ann.score_edges(edges, regulon, genes, use_sign=False)
            signed = ann.score_edges(edges, regulon, genes, use_sign=True)
            enr = ann.enrichment(edges, regulon, genes, n_draws=args.draws)
            entry["sources"][name] = {
                "n_edges": len(edges), "unsigned": unsigned, "signed": signed,
                "enrichment": enr, "delta_vs_linear": deltas.get(name),
            }
            p = unsigned["precision"]
            print(f"  {name:<16}{len(edges):>4} edges  hits {unsigned['n_hit']:>3}  "
                  f"precision {p:.4f}  chance {enr['chance_precision']:.4f}  "
                  f"ratio {enr['ratio']:.2f}  p {enr['p']:.4f}  "
                  f"signed hits {signed['n_hit']:>3}")
        out["modules"][module] = entry

    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nwrote {outpath}")


if __name__ == "__main__":
    main()
