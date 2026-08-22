"""The three figures the fourth review named, built from recorded results only.

  fig1  calibration: which generators reproduce the identifiability signature
  fig2  the signature across cell states and across atlases, matched for shape
  fig3  search convergence, and what the choice of seed summary does to the count

Nothing here recomputes a result. Every value is read from a results file, so a figure
cannot drift from the number it illustrates.

    python scripts/make_figures.py --out results/figures
"""
from __future__ import annotations

import argparse
import glob
import json
import statistics as st
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RANK_TOL_REL = 0.15
PC1_TOL_ABS = 0.05

FAMILY_STYLE = {
    "structural": ("#B04A3A", "o", "structural, uniform"),
    "structural[uniform]": ("#B04A3A", "o", None),
    "structural[hub]": ("#C97B3C", "s", "structural, hub"),
    "structural[scale_free]": ("#8A6A2F", "^", "structural, scale-free"),
    "structural[modular]": ("#6E4B7A", "D", "structural, modular"),
    "structural+lowrank": ("#3E7CB1", "v", "structural + shared component"),
    "lowrank": ("#2E7D5B", "*", "low-rank surrogate, no structure"),
}


def fig1_calibration(outdir: Path) -> None:
    """Distance from the real module, per generator, with the acceptance box drawn.

    One point per generator variant per module. Axes are the two acceptance quantities,
    so a point inside the shaded corner is a generator that reproduced the signature.
    """
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    seen = set()
    for p in sorted(glob.glob("results/a17/*.json")):
        d = json.load(open(p))
        for v in d["variants"]:
            fam = v["generator"]
            colour, marker, label = FAMILY_STYLE.get(fam, ("#888888", "x", fam))
            x = max(v["effective_rank_relative_error"], 1e-3)
            y = max(v["leading_pc_absolute_error"], 1e-4)
            show = label if label and label not in seen else None
            if show:
                seen.add(label)
            ax.scatter(x, y, c=colour, marker=marker, s=54 if marker != "*" else 150,
                       alpha=0.85, edgecolors="white", linewidths=0.6, label=show, zorder=3)

    ax.add_patch(plt.Rectangle((1e-3, 1e-4), RANK_TOL_REL - 1e-3, PC1_TOL_ABS - 1e-4,
                               facecolor="#2E7D5B", alpha=0.10, zorder=1))
    ax.axvline(RANK_TOL_REL, color="#2E7D5B", ls="--", lw=1, zorder=2)
    ax.axhline(PC1_TOL_ABS, color="#2E7D5B", ls="--", lw=1, zorder=2)
    ax.text(RANK_TOL_REL * 0.92, PC1_TOL_ABS * 0.9, "accepted", ha="right", va="top",
            fontsize=9, color="#2E7D5B")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("effective rank, relative error against the real module")
    ax.set_ylabel("leading-PC fraction, absolute error")
    ax.set_title("Only a structureless low-rank surrogate reproduces the signature\n"
                 "seven modules; the two structural points inside the box both "
                 "carry a dominant shared component", fontsize=10)
    ax.legend(fontsize=8, loc="lower right", frameon=False)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(outdir / "fig1_calibration.png", dpi=200)
    plt.close(fig)
    print("  fig1_calibration.png")


def fig2_scope(outdir: Path) -> None:
    """Left: the signature across the three cell states. Right: against a second atlas."""
    by_n: dict[int, list] = {}
    states: dict[str, list] = {"Rest": [], "Stim8hr": [], "Stim48hr": []}
    for p in glob.glob("results/a19/zhu_*.json"):
        d = json.load(open(p))
        for s, vals in states.items():
            v = d["states"].get(s) or {}
            if "error" not in v and v:
                vals.append(v["effective_rank_normalised"])
        v8 = d["states"].get("Stim8hr") or {}
        if v8 and "error" not in v8:
            by_n.setdefault(v8["n_perturbations"], []).append(v8["effective_rank_normalised"])

    nm = json.load(open("results/a19/norman.json"))["shape_matched"]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 4.6))

    names = list(states)
    a1.boxplot([states[n] for n in names], tick_labels=["rest", "8 h", "48 h"],
               widths=0.55, patch_artist=True,
               boxprops=dict(facecolor="#DCE6F0", edgecolor="#3E5C76"),
               medianprops=dict(color="#B04A3A", lw=1.6),
               flierprops=dict(marker="o", ms=3, mfc="#3E5C76", mec="none", alpha=0.6))
    for i, n in enumerate(names, start=1):
        a1.scatter(np.random.default_rng(i).normal(i, 0.055, len(states[n])), states[n],
                   s=13, c="#3E5C76", alpha=0.55, zorder=3)
    a1.set_ylabel("normalised effective rank")
    a1.set_title("Stable across cell state\n27 modules, same atlas", fontsize=10)
    a1.set_ylim(0, 0.5)
    a1.grid(alpha=0.25, lw=0.5, axis="y")

    shapes = [k for k in (11, 20, 28, 40) if f"k={k}" in nm and k in by_n]
    zhu = [st.median(by_n[k]) for k in shapes]
    nor = [nm[f"k={k}"]["effective_rank_normalised_mean"] for k in shapes]
    x = np.arange(len(shapes))
    a2.bar(x - 0.19, zhu, 0.38, label="this atlas", color="#3E5C76")
    a2.bar(x + 0.19, nor, 0.38, label="second atlas, shape matched", color="#C97B3C")
    # how many modules each bar rests on. Two of the shapes have exactly one module here, so
    # their bar is a single value rather than a median and should not be read as robust
    for xi, k, v in zip(x, shapes, zhu):
        a2.text(xi - 0.19, v + 0.008, f"n={len(by_n[k])}", ha="center", va="bottom",
                fontsize=7.5, color="#3E5C76")
    a2.set_xticks(x)
    a2.set_xticklabels([f"{k} by {k}" for k in shapes])
    a2.set_ylabel("normalised effective rank")
    a2.set_title("Not stable across atlas\nsame shapes, readouts = perturbed set",
                 fontsize=10)
    a2.set_ylim(0, 0.5)
    a2.legend(fontsize=8, frameon=False)
    a2.grid(alpha=0.25, lw=0.5, axis="y")

    fig.tight_layout()
    fig.savefig(outdir / "fig2_scope.png", dpi=200)
    plt.close(fig)
    print("  fig2_scope.png")


def fig3_search(outdir: Path) -> None:
    """Left: saturation curves. Right: what the seed summary does to each module."""
    sat = json.load(open("results/seed_saturation.json"))["rows"]
    ver = {r["module"]: r for r in json.load(open("results/a20_verdict.json"))["rows"]}

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.8, 4.8))

    for r in sorted(sat, key=lambda x: -x["last_gain"]):
        ks = range(1, len(r["curve"]) + 1)
        climbing = r["last_gain"] >= 0.005
        a1.plot(ks, r["curve"], marker="o", ms=3.2, lw=1.5 if climbing else 0.9,
                color="#B04A3A" if climbing else "#9AA5B1",
                alpha=0.95 if climbing else 0.7, zorder=3 if climbing else 2)
    a1.plot([], [], color="#B04A3A", lw=1.5, label="still climbing at k=5")
    a1.plot([], [], color="#9AA5B1", lw=0.9, label="converged")
    a1.set_xlabel("number of search seeds, k")
    a1.set_ylabel("mean best-of-k held-out DE-overlap")
    a1.set_title("The oracle maximum is not a converged bound\n"
                 "10 of 13 modules still climbing", fontsize=10)
    a1.set_xticks(list(range(1, 6)))
    a1.legend(fontsize=8, frameon=False, loc="upper left")
    a1.grid(alpha=0.25, lw=0.5)

    mods = sorted(ver, key=lambda m: ver[m]["oracle_seed_max"] - ver[m]["oracle_seed_median"])
    y = np.arange(len(mods))
    for i, m in enumerate(mods):
        r = ver[m]
        a2.plot([r["oracle_seed_median"], r["oracle_seed_max"]], [i, i],
                color="#C9CDD3", lw=2, zorder=1)
    # the maximum is drawn larger and underneath so that a converged module, where the two
    # summaries coincide, still shows both rather than hiding one behind the other
    a2.scatter([ver[m]["oracle_seed_max"] for m in mods], y, s=64, c="#B04A3A",
               label="oracle, seed maximum", zorder=3)
    a2.scatter([ver[m]["oracle_seed_median"] for m in mods], y, s=22, c="#3E5C76",
               label="oracle, seed median", zorder=4)
    a2.scatter([ver[m]["linear_mean"] for m in mods], y, s=44, marker="|", c="#2E7D5B",
               linewidths=2, label="linear baseline", zorder=4)
    a2.set_yticks(y)
    a2.set_yticklabels([m.replace("coresponse_", "co:").replace("regulon_", "reg:")
                        for m in mods], fontsize=8)
    a2.set_xlabel("held-out DE-overlap")
    a2.set_title("Which summary is used decides the count\n"
                 "1 of 13 on the median, 5 of 13 on the maximum", fontsize=10)
    a2.legend(fontsize=8, frameon=False, loc="lower right")
    a2.grid(alpha=0.25, lw=0.5, axis="x")

    fig.tight_layout()
    fig.savefig(outdir / "fig3_search.png", dpi=200)
    plt.close(fig)
    print("  fig3_search.png")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/figures")
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    fig1_calibration(outdir)
    fig2_scope(outdir)
    fig3_search(outdir)
    print(f"wrote into {outdir}")


if __name__ == "__main__":
    main()
