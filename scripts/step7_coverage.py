"""Coverage: how often does a sparse structure actually make a prediction?

A structural prediction is the difference between the clamped and unclamped fixed points,
so a perturbation of a gene with no outgoing edges moves nothing and the model predicts
identically zero for every gene. Those folds are scored anyway, at chance, which means a
structure's headline score mixes folds where it says something with folds where it says
nothing, while the linear baseline predicts on every fold.

This reports coverage per source as folds predicted over folds scoreable, and repeats the
paired comparison against linear on the covered subset alone. If a structure still loses
where it applies, the architectural forfeit is not what is holding it back.

Pure re-analysis of existing comparator results; nothing is re-run.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

from mmc.eval import compare


def regulators(spec: dict) -> set[str]:
    out = set()
    for e in spec.get("edges", []):
        if isinstance(e, dict):
            r = e.get("regulator") or e.get("source")
        else:
            r = e[0]
        if r:
            out.add(r)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    for pattern in args.results:
        for f in sorted(glob.glob(pattern)):
            d = json.load(open(f))
            if not (d.get("table") and d.get("random_null")):
                continue
            perts = d["module_summary"]["perts"]
            pf = d.get("per_fold_scores") or {}
            if "linear" not in pf:
                continue
            for src, spec in (d.get("specs") or {}).items():
                if src not in pf or not isinstance(spec, dict):
                    continue
                regs = regulators(spec)
                idx = [i for i in range(len(perts))
                       if pf[src][i] is not None and pf["linear"][i] is not None]
                cov = [i for i in idx if perts[i] in regs]
                if not idx:
                    continue
                a_all = np.array([pf[src][i] for i in idx], float)
                b_all = np.array([pf["linear"][i] for i in idx], float)
                rec = {"module": d["module"], "source": src,
                       "n_edges": len(spec.get("edges", [])),
                       "n_regulators": len(regs),
                       "folds_scoreable": len(idx), "folds_covered": len(cov),
                       "coverage": len(cov) / len(idx),
                       "all_folds": compare.paired_delta(a_all, b_all),
                       "mean_on_uncovered": float(np.mean(
                           [pf[src][i] for i in idx if i not in set(cov)]))
                       if len(cov) < len(idx) else None}
                if len(cov) >= 3:
                    a_c = np.array([pf[src][i] for i in cov], float)
                    b_c = np.array([pf["linear"][i] for i in cov], float)
                    rec["covered_folds"] = compare.paired_delta(a_c, b_c)
                    rec["source_mean_covered"] = float(a_c.mean())
                    rec["linear_mean_covered"] = float(b_c.mean())
                else:
                    rec["covered_folds"] = None
                rows.append(rec)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump({"rows": rows}, open(args.out, "w"), indent=2, default=float)

    ora = [r for r in rows if r["source"] == "oracle"]
    cov = np.array([r["coverage"] for r in ora], float)
    print(f"{len(ora)} modules with an oracle structure\n")
    print(f"coverage, folds where the perturbed gene has outgoing edges:")
    print(f"  median {np.median(cov):.2f}   min {cov.min():.2f}   max {cov.max():.2f}"
          f"   below half: {int((cov < 0.5).sum())} of {len(cov)}\n")
    scored = [r for r in ora if r["covered_folds"]]
    win_all = sum(1 for r in ora if r["all_folds"]["advantage"])
    win_cov = sum(1 for r in scored if r["covered_folds"]["advantage"])
    print(f"oracle clears linear on all scoreable folds:   {win_all} of {len(ora)}")
    print(f"oracle clears linear on covered folds alone:   {win_cov} of {len(scored)}"
          f"   (modules with at least 3 covered folds)\n")
    print(f"{'module':<24}{'cov':>6}{'edges':>7}{'delta all':>12}{'delta covered':>16}")
    for r in sorted(ora, key=lambda r: r["coverage"]):
        c = r["covered_folds"]
        print(f"  {r['module']:<22}{r['coverage']:>6.2f}{r['n_edges']:>7}"
              f"{r['all_folds']['delta']:>+12.4f}"
              + (f"{c['delta']:>+16.4f}" if c else f"{'too few':>16}"))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
