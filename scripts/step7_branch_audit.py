"""Does the ceiling's advantage over linear survive where linear is worth beating?

Family A counts modules where the oracle ceiling beats the linear baseline. Read alone
that count says structure wins somewhere, which sits badly beside the finding that no
structure source beats linear. The resolution is that the linear arm is at or below its
own random null on a large fraction of modules, and beating a baseline that carries no
signal is not evidence of anything.

This cross-tabulates the two facts over every module: whether the oracle clears linear
after the Benjamini-Hochberg correction of section 1.3, against whether linear clears its
own random null. Only the cell where both hold can bear on whether Branch A generalises.

Pure re-analysis of existing comparator results; nothing is re-run.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path


def null_reference(d: dict) -> tuple[float, float, str]:
    """The null this module's scores are judged against: mean and 95th percentile.

    Swept nulls carry one band per edge count; the band matching the oracle's structure
    is the comparable one, and the widest band is used when no match is recorded.
    """
    bands = d.get("random_null_bands") or {}
    spec = (d.get("specs") or {}).get("oracle") or {}
    n_edges = None
    if isinstance(spec, dict):
        n_edges = spec.get("n_edges")
    if not n_edges:
        tab = {r["source"]: r for r in d.get("table", [])}
        n_edges = (tab.get("oracle") or {}).get("n_edges")
    if bands:
        key = str(n_edges) if str(n_edges) in bands else sorted(bands)[-1]
        b = bands[key]
        return float(b["mean"]), float(b.get("p95", b["mean"])), f"band {key}"
    rn = d.get("random_null") or {}
    return float(rn.get("mean", float("nan"))), float(rn.get("p95", float("nan"))), "single"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True, help="step2_regime_fit.py output")
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fit = json.load(open(args.fit))
    by_name = {r["module"]: r for r in fit["rows"]}

    raw = {}
    for pattern in args.results:
        for f in sorted(glob.glob(pattern)):
            d = json.load(open(f))
            if d.get("table") and d.get("random_null"):
                raw[d["module"]] = d

    rows = []
    for name, r in by_name.items():
        d = raw.get(name)
        if d is None:
            continue
        null_mean, null_p95, kind = null_reference(d)
        linear = r["linear"]
        rows.append({
            "module": name, "source": r["source"], "n_folds": r["n_folds"],
            "linear": linear, "oracle": r["oracle"], "null_mean": null_mean,
            "null_p95": null_p95, "null_kind": kind,
            "oracle_clears_linear": bool(r.get("bh_reject")),
            "oracle_nominal": bool(r.get("advantage")),
            "linear_beats_null": bool(linear > null_p95),
            "linear_above_null_mean": bool(linear > null_mean),
        })

    cells = {}
    for a in (True, False):
        for b in (True, False):
            cells[f"oracle_clears_linear={a},linear_beats_null={b}"] = [
                x["module"] for x in rows
                if x["oracle_clears_linear"] == a and x["linear_beats_null"] == b]

    out = {"n_modules": len(rows), "rows": rows, "cells": cells,
           "linear_beats_null": sum(1 for x in rows if x["linear_beats_null"]),
           "oracle_clears_linear": sum(1 for x in rows if x["oracle_clears_linear"])}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=2)

    n = len(rows)
    tt = len(cells["oracle_clears_linear=True,linear_beats_null=True"])
    tf = len(cells["oracle_clears_linear=True,linear_beats_null=False"])
    ft = len(cells["oracle_clears_linear=False,linear_beats_null=True"])
    ff = len(cells["oracle_clears_linear=False,linear_beats_null=False"])
    print(f"{n} modules. Linear is judged to carry signal when it exceeds the 95th")
    print("percentile of its module's own random null.\n")
    print(f"{'':<34}{'linear carries signal':>22}{'linear is noise':>18}")
    print(f"{'oracle clears linear (post-BH)':<34}{tt:>22}{tf:>18}")
    print(f"{'oracle does not clear linear':<34}{ft:>22}{ff:>18}")
    print(f"\nthe only cell that can refute Branch A is the top left, n = {tt}")
    for label, key in (("oracle clears linear, linear carries signal",
                        "oracle_clears_linear=True,linear_beats_null=True"),
                       ("oracle clears linear, linear is noise (vacuous)",
                        "oracle_clears_linear=True,linear_beats_null=False")):
        ms = cells[key]
        print(f"\n{label} ({len(ms)}):")
        for m in ms:
            x = next(r for r in rows if r["module"] == m)
            print(f"  {m:<24} linear {x['linear']:.4f}  null p95 {x['null_p95']:.4f}  "
                  f"oracle {x['oracle']:.4f}  {x['n_folds']} folds")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
