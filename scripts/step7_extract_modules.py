"""Step 7: extract and screen candidate modules for the scale-up (PREREG_v4 section 8).

Builds candidates from two sources, screens each against the power precondition, and
writes **every** candidate with its verdict. Section 8 requires that modules failing the
screen are reported rather than dropped, so that the set entering the regime map cannot
have been curated toward the conclusion.

The screen's fold-count floor is the lesson from Th2_GATA3, which ran with two scoreable
folds and produced intervals spanning the entire metric. A module that cannot support the
paired statistic should be identified before it is run, not after.

    python scripts/step7_extract_modules.py --collectri /data/collectri.tsv \\
        --condition Stim8hr --out /work/results/step7_modules.json

Requires MMC_ZHU_STORE. Extraction only; running the comparator on the passing modules
is a separate job.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mmc.data import module_generator as mg
from mmc.data.module_data import build_module_data
from mmc.data import module_extract
from mmc.eval import annotation as ann
from mmc.shared import store


def characterise_candidate(cand: mg.Candidate, condition: str) -> mg.Candidate:
    """Assemble the candidate's response matrix and record its shape."""
    name = f"cand_{cand.name}"
    module_extract.register_module(name, cand.genes, cand.genes)
    try:
        mod = build_module_data(name, condition)
    except ValueError:
        cand.n_perts = 0
        return mg.screen(cand)
    stats = mg.characterise(cand.genes, mod.observed, mod.de_mask)
    cand.n_perts = stats["n_perts"]
    cand.n_de_entries = stats["n_de_entries"]
    cand.n_scoreable_folds = stats["n_scoreable_folds"]
    return mg.screen(cand)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--collectri", required=True)
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--universe", type=int, default=400,
                    help="most-perturbed genes to consider for co-response clustering")
    ap.add_argument("--regulon-modules", type=int, default=40)
    ap.add_argument("--coresponse-modules", type=int, default=20)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    measured = store.measured_genes()
    con = store._con()
    perturbed = {str(r[0]) for r in con.execute(
        "SELECT DISTINCT perturbation FROM zhu_pert WHERE condition = ?",
        [args.condition]).fetchall()}
    print(f"atlas: {len(measured)} measured genes, {len(perturbed)} perturbed in "
          f"{args.condition}")

    regulon = ann.load_collectri(args.collectri)
    print(f"CollecTRI: {len(regulon.edges)} interactions")

    cands = mg.regulon_candidates(regulon, measured, perturbed,
                                  max_modules=args.regulon_modules)
    print(f"regulon candidates: {len(cands)}")

    # co-response clustering over the most broadly-acting perturbed genes, which are the
    # ones a module can actually be exercised on
    top = con.execute(
        "SELECT perturbation, MAX(n_downstream) AS breadth FROM zhu_pert "
        "WHERE condition = ? GROUP BY perturbation ORDER BY breadth DESC LIMIT ?",
        [args.condition, args.universe]).df()
    universe = [g for g in top["perturbation"].astype(str).tolist() if g in measured]
    print(f"co-response universe: {len(universe)} genes")

    if len(universe) >= 20:
        uni_name = "cand_universe"
        module_extract.register_module(uni_name, universe, universe)
        try:
            uni = build_module_data(uni_name, args.condition)
            cands += mg.coresponse_candidates(list(uni.genes), uni.observed,
                                              size=20, n_modules=args.coresponse_modules)
            print(f"co-response candidates: {args.coresponse_modules} requested")
        except ValueError as e:
            print(f"co-response clustering skipped: {e}")

    print(f"\nscreening {len(cands)} candidates")
    screened = []
    for i, c in enumerate(cands, 1):
        screened.append(characterise_candidate(c, args.condition))
        mark = "PASS" if c.passed else "fail"
        print(f"  [{i:>3}/{len(cands)}] {c.name:<28} {len(c.genes):>3} genes  "
              f"{c.n_perts:>3} perts  {c.n_de_entries:>4} DE  "
              f"{c.n_scoreable_folds:>3} folds  {mark}"
              + ("" if c.passed else f"  ({'; '.join(c.reasons)})"))

    summary = mg.summarise(screened)
    out = {"condition": args.condition, "summary": summary,
           "candidates": [c.as_dict() for c in screened]}
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2, default=float)

    print(f"\nextracted {summary['n_extracted']}, passed {summary['n_passed']}, "
          f"failed {summary['n_failed']}")
    print(f"by source: {summary['by_source']}")
    print(f"passed by source: {summary['passed_by_source']}")
    if summary["failure_reasons"]:
        print("failure reasons:")
        for r, n in sorted(summary["failure_reasons"].items(), key=lambda kv: -kv[1]):
            print(f"  {n:>3}  {r}")
    passed = [c for c in screened if c.passed]
    if passed:
        folds = np.array([c.n_scoreable_folds for c in passed])
        print(f"passing modules carry {folds.min()} to {folds.max()} scoreable folds "
              f"(median {int(np.median(folds))})")
    print(f"wrote {outpath}")


if __name__ == "__main__":
    main()
