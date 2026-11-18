"""Pool the offset arms across modules: does the offset lift any structure, or only a good one?

One reading of a structural model losing to a linear map is that the model class is the
constraint, and that the structure it was given would predict if the class could express the
shared response. Adding a fitted per-gene offset tests that reading, and on a single searched
structure the offset does move the arm past the linear map.

The test only means something against structures of known quality, so the offset is applied to
three per module, all scored on the real responses under one protocol:

  leaky      the structure from the twenty-seed search, selected on the folds it is scored on
  nested     the structure from the nested search, selected on inner folds only
  permuted   the structure from the search on data with the perturbation-to-response pairing
             destroyed, scored here against the real responses, so by construction it carries
             no information about them

The permuted arm is the one that decides it. If an offset lifts a structure that cannot contain
real information by as much as it lifts a searched one, the lift is a property of the offset
rather than evidence that a missing structure was holding the grammar back.

Reads the per-module files `run_offset_arm.py` writes and refits nothing.

    python scripts/summarize_offset_arms.py \
        --dir results/comparator/offset \
        --out results/comparator/offset_summary.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

import numpy as np


def score(table, name):
    for r in table:
        if r["source"] == name:
            return r["de_overlap_mean"]
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/comparator/offset")
    ap.add_argument("--out", default="results/comparator/offset_summary.json")
    args = ap.parse_args()

    rows = []
    for f in sorted(glob.glob(f"{args.dir}/real_*.json")):
        d = json.load(open(f))
        m = d["module"]
        t = d["table"]
        pf = os.path.join(args.dir, f"permuted_{m}.json")
        pt = json.load(open(pf))["table"] if os.path.exists(pf) else []
        row = {
            "module": m,
            "linear": score(t, "linear"),
            "leaky": score(t, "oracle"),
            "leaky_offset": score(t, "oracle+offset"),
            "nested": score(t, "oracle_nested"),
            "nested_offset": score(t, "oracle_nested+offset"),
            "permuted": score(pt, "oracle"),
            "permuted_offset": score(pt, "oracle+offset"),
        }
        for arm in ("leaky", "nested", "permuted"):
            base, off = row[arm], row[f"{arm}_offset"]
            row[f"{arm}_lift"] = (off - base) if (base is not None and off is not None) else None
            row[f"{arm}_offset_beats_linear"] = bool(
                off is not None and row["linear"] is not None and off > row["linear"])
        rows.append(row)

    print(f"{'module':<24}{'linear':>8} | {'leaky':>16} | {'nested':>16} | {'permuted':>16}")
    print(f"{'':<24}{'':>8} | {'base -> offset':>16} | {'base -> offset':>16} | "
          f"{'base -> offset':>16}")
    for r in rows:
        def cell(a, r=r):
            b, o = r[a], r[f"{a}_offset"]
            if b is None or o is None:
                return f"{'n/a':>16}"
            return f"{b:>7.4f}->{o:<8.4f}"
        print(f"{r['module']:<24}{r['linear']:>8.4f} | {cell('leaky')} | {cell('nested')} | "
              f"{cell('permuted')}")

    print(f"\n{'arm':<12}{'lifts past linear':>20}{'median lift':>14}")
    summary = {}
    for arm in ("leaky", "nested", "permuted"):
        beats = sum(r[f"{arm}_offset_beats_linear"] for r in rows)
        lifts = [r[f"{arm}_lift"] for r in rows if r[f"{arm}_lift"] is not None]
        summary[arm] = {"offset_beats_linear": beats, "n": len(rows),
                        "median_lift": float(np.median(lifts)) if lifts else None,
                        "n_lifts": len(lifts)}
        med = f"{np.median(lifts):+.4f}" if lifts else "n/a"
        print(f"{arm:<12}{f'{beats} of {len(rows)}':>20}{med:>14}")

    print("\nThe permuted structure cannot contain information about the responses it is scored\n"
          "against, so whatever the offset buys it is what an offset buys any structure of that\n"
          "size.")

    out = {"n_modules": len(rows), "summary": summary, "rows": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
