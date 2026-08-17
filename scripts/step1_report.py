"""Assemble the Step 1 comparator tables and figure from the per-module result files.

Reads the JSON written by `step1_comparator.py` for each module and produces the
combined table, the Benjamini-Hochberg decision within Family A, the structure
agreement matrix, and Figure 1.

    python scripts/step1_report.py results/step1_*.json --out paper/step1

One thing this script insists on recording, because the comparator is easy to
misread without it: where each structure source's edges came from relative to the
held-out folds. The harness refits parameters inside every fold, but the structure
itself is fixed across folds for every source. So a source that derived its edges
from the full response matrix has seen the held-out perturbations at selection time,
even though it never sees them at fit time. That is true of the oracle by design, and
it is equally true of the mean-difference and GRNBoost2 arms, which read the same
matrix. The textbook arm is the only one whose edges were fixed without reference to
this data at all. The nested oracle is the one arm with an outer split that selection
never touched. Reporting a ranking without this column would invite the reader to
treat all the rows as comparable in a way they are not.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

from mmc.eval import compare

# How each source's edge set was obtained relative to the evaluated folds.
PROVENANCE = {
    "textbook": ("prior", "fixed from published biology before this data was read"),
    "claude": ("module", "proposed from the module's training data at loop time"),
    "mean_difference": ("full matrix", "thresholded from the full response matrix"),
    "grnboost2": ("full matrix", "regressed on the full response matrix"),
    "oracle": ("held-out", "selected on the held-out score; ceiling estimator"),
    "oracle_nested": ("inner split", "selected on an inner split, scored on an untouched outer split"),
    "linear": ("not a structure", "the bar"),
    "mean": ("not a structure", "mean of training perturbations"),
    "zero": ("not a structure", "predicts no change"),
}
STRUCTURE_SOURCES = ("textbook", "claude", "mean_difference", "grnboost2", "oracle")


def load(paths: list[str]) -> dict[str, dict]:
    out = {}
    for p in sorted(paths):
        with open(p) as f:
            d = json.load(f)
        out[d["module"]] = d
    return out


def family_a(results: dict[str, dict], source: str) -> list[dict]:
    """The module-level test that `source` beats linear, across modules (PREREG 1.3)."""
    rows = []
    for module, d in results.items():
        per_fold = d.get("per_fold_scores", {})
        if source not in per_fold or "linear" not in per_fold:
            continue
        a = np.array([np.nan if v is None else v for v in per_fold[source]], float)
        b = np.array([np.nan if v is None else v for v in per_fold["linear"]], float)
        delta = compare.paired_delta(a, b)
        rows.append({"module": module, "source": source,
                     "delta": round(delta["delta"], 4),
                     "lo": round(delta["lo"], 4), "hi": round(delta["hi"], 4),
                     "n_folds": delta["n_folds"],
                     "p": round(compare.permutation_p(a, b), 4),
                     "advantage": delta["advantage"]})
    if rows:
        flags = compare.benjamini_hochberg([r["p"] for r in rows], q=0.05)
        for r, f in zip(rows, flags):
            r["bh_reject"] = bool(f)
    return rows


def branch_read(results: dict[str, dict], primary: str = "Cytokine_production") -> dict:
    """The PREREG_v4 section 2.4 branch, read on the primary module."""
    d = results.get(primary)
    if not d:
        return {"error": f"{primary} result not present"}
    pf = d.get("per_fold_scores", {})
    if "linear" not in pf:
        return {"error": "no linear arm"}
    lin = np.array([np.nan if v is None else v for v in pf["linear"]], float)

    def delta(name):
        if name not in pf:
            return None
        arr = np.array([np.nan if v is None else v for v in pf[name]], float)
        return compare.paired_delta(arr, lin)

    oracle, claude = delta("oracle"), delta("claude")
    out = {"module": primary,
           "oracle_vs_linear": oracle,
           "claude_vs_linear": claude,
           "nested_outer": d.get("oracle_nested", {}).get("outer_de_overlap"),
           "ceiling_spread": d.get("oracle_ceiling_spread")}
    if oracle is None:
        out["branch"] = "undetermined: no oracle arm"
    elif claude is not None and claude["advantage"]:
        out["branch"] = "C"
        out["reading"] = ("the proposal arm clears the linear baseline; the pre-registration "
                          "requires this be honoured and the earlier module-level negatives "
                          "re-read as underpowered")
    elif oracle["advantage"]:
        out["branch"] = "B"
        out["reading"] = ("the data supports a predictive structure that the proposer does not "
                          "find; the failure localises to proposal and search")
    else:
        out["branch"] = "A"
        out["reading"] = ("no structure in this grammar, including one selected on the held-out "
                          "answer, beats a linear map; the ceiling is set by the data")
    return out


def figure(results: dict[str, dict], out_prefix: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    modules = list(results)
    fig, axes = plt.subplots(1, len(modules), figsize=(4.2 * len(modules), 4.6), squeeze=False)
    for ax, module in zip(axes[0], modules):
        d = results[module]
        table = d.get("table", [])
        names = [r["source"] for r in table]
        vals = [r["de_overlap_mean"] for r in table]
        los = [r["de_overlap_mean"] - r["de_overlap_lo"] for r in table]
        his = [r["de_overlap_hi"] - r["de_overlap_mean"] for r in table]
        colors = ["#b03030" if n == "oracle" else
                  "#30609a" if n in ("linear", "mean", "zero") else "#7a7a7a"
                  for n in names]
        y = np.arange(len(names))
        ax.barh(y, vals, xerr=[los, his], color=colors, height=0.65,
                error_kw={"elinewidth": 1, "capsize": 2})
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.invert_yaxis()
        null = d.get("random_null", {})
        if null.get("mean") is not None:
            ax.axvline(null["mean"], color="#444", ls="--", lw=1)
            ax.axvline(null.get("p95", null["mean"]), color="#444", ls=":", lw=1)
        summ = d.get("module_summary", {})
        ax.set_title(f"{module}\n{summ.get('n_perts', '?')} perturbations, "
                     f"{summ.get('n_de_entries', '?')} DE entries", fontsize=9)
        ax.set_xlabel("held-out DE-overlap")
        ax.grid(axis="x", alpha=0.25)
    fig.suptitle("Structure sources against the linear baseline. Dashed line: random-structure "
                 "null mean; dotted: its 95th percentile.", fontsize=9, y=0.02)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(f"{out_prefix}_fig1.png", dpi=200)
    print(f"wrote {out_prefix}_fig1.png")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--out", default="paper/step1")
    ap.add_argument("--primary", default="Cytokine_production")
    args = ap.parse_args()

    paths = [p for pattern in args.paths for p in glob.glob(pattern)]
    results = load(paths)
    if not results:
        raise SystemExit("no result files matched")
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    report = {
        "modules": {m: d.get("module_summary") for m, d in results.items()},
        "protocol": next(iter(results.values())).get("protocol"),
        "provenance": {k: {"selected_from": v[0], "note": v[1]} for k, v in PROVENANCE.items()},
        "tables": {m: d.get("table") for m, d in results.items()},
        "diagnostics": {m: d.get("diagnostics") for m, d in results.items()},
        "random_null": {m: {k: v for k, v in (d.get("random_null") or {}).items()
                            if k != "values"} for m, d in results.items()},
        "structure_agreement": {m: d.get("structure_agreement") for m, d in results.items()},
        "family_a": {s: family_a(results, s) for s in STRUCTURE_SOURCES},
        "branch": branch_read(results, args.primary),
        "algorithmic_scope": next((d.get("algorithmic_scope") for d in results.values()
                                   if d.get("algorithmic_scope")), None),
    }
    with open(f"{out_prefix}_report.json", "w") as f:
        json.dump(report, f, indent=2, default=float)
    print(f"wrote {out_prefix}_report.json")

    for module, d in results.items():
        print(f"\n=== {module} ===")
        s = d.get("module_summary", {})
        print(f"{s.get('n_genes')} genes, {s.get('n_perts')} perturbations, "
              f"{s.get('n_de_entries')} DE entries, {s.get('n_folds_with_de')} scoreable folds")
        print(f"  {'source':<18}{'edges':>6}  {'DE-overlap [95% CI]':<26}"
              f"{'delta vs linear':<26}{'null pct':>9}  {'selected from':<14}")
        for r in d.get("table", []):
            prov = PROVENANCE.get(r["source"], ("", ""))[0]
            pct = None if prov == "not a structure" else r.get("random_null_percentile")
            print(f"  {r['source']:<18}{r['n_edges']:>6}  "
                  f"{r['de_overlap_mean']:.4f} [{r['de_overlap_lo']:.4f}, {r['de_overlap_hi']:.4f}]  "
                  f"{r['delta_vs_linear']:+.4f} [{r['delta_lo']:+.4f}, {r['delta_hi']:+.4f}]  "
                  f"{'' if pct is None else f'{pct:>7.1f}%'}  {prov:<14}")

    b = report["branch"]
    print(f"\n=== branch (read on {b.get('module')}) ===")
    if "branch" in b:
        print(f"  Branch {b['branch']}: {b.get('reading', '')}")
        o = b.get("oracle_vs_linear")
        if o:
            print(f"  oracle minus linear: {o['delta']:+.4f} [{o['lo']:+.4f}, {o['hi']:+.4f}] "
                  f"over {o['n_folds']} folds")
        sp = b.get("ceiling_spread")
        if sp:
            print(f"  ceiling across seeds: {sp['min']:.4f} to {sp['max']:.4f} (sd {sp['sd']:.4f})")
        ns = b.get("nested_outer")
        if ns:
            print(f"  nested oracle on the untouched outer split: {ns['mean']:.4f} "
                  f"[{ns['lo']:.4f}, {ns['hi']:.4f}]")
    else:
        print(" ", b.get("error"))

    try:
        figure(results, out_prefix)
    except Exception as e:                                    # noqa: BLE001
        print(f"figure skipped: {e}")


if __name__ == "__main__":
    main()
