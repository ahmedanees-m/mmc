"""Render the scaled engineer-behavior figure (two panels).

Panel A: the held-out validated rate of novel hypotheses and the module-condition
beats-linear rate, each as a point estimate with a Wilson score confidence interval. Panel B:
the reasoning-versus-search comparison, the proposed structure against a random structure of
equal edge count on in-sample fit, per module.

Reads engineer_behavior_scaled.json.
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

GREEN, GREY, RED = "#2E7D32", "#9E9E9E", "#C62828"


def wilson(k, n):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    z = 1.96
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else "engineer_behavior_scaled.json"
    out = sys.argv[2] if len(sys.argv) > 2 else "engineer_behavior_scaled.png"
    with open(path) as f:
        d = json.load(f)

    ne = d["novel_edges"]
    distinct = {e["edge"] for e in ne}
    dv = {e["edge"] for e in ne if e["in_beats_linear_module"]}
    beat, nmod = d["summary"]["modules_beating_linear"]

    rows = [
        (f"novel hypotheses in a model\nthat beats linear ({len(dv)} of {len(distinct)})",
         wilson(len(dv), len(distinct))),
        (f"module-conditions that\nbeat linear held-out ({beat} of {nmod})",
         wilson(beat, nmod)),
    ]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 4.6))

    y = range(len(rows))
    for i, (_lab, (p, lo, hi)) in enumerate(rows):
        axA.plot([lo, hi], [i, i], color=GREY, lw=3, zorder=1)
        axA.scatter([p], [i], color=RED, zorder=2, s=60)
        axA.text(hi + 0.01, i, f"95% CI [{lo*100:.0f}%, {hi*100:.0f}%]", va="center", fontsize=9)
    axA.set_yticks(list(y))
    axA.set_yticklabels([r[0] for r in rows], fontsize=9)
    axA.set_xlim(-0.02, 0.5)
    axA.set_xlabel("held-out rate (point estimate and Wilson 95% CI)")
    axA.set_title("A  No held-out advantage over a linear baseline", fontsize=11)

    mods = defaultdict(lambda: [[], []])
    for m in d["modules"]:
        if m.get("claude_insample_pearson") is not None:
            mods[m["module"]][0].append(m["claude_insample_pearson"])
            mods[m["module"]][1].append(m["random_insample_pearson"])
    names = list(mods.keys())
    claude = [sum(mods[k][0]) / len(mods[k][0]) for k in names]
    rand = [sum(mods[k][1]) / len(mods[k][1]) for k in names]
    x = range(len(names))
    axB.bar([i - 0.2 for i in x], claude, width=0.4, color=GREEN, label="reasoning (Claude)")
    axB.bar([i + 0.2 for i in x], rand, width=0.4, color=GREY, label="random (equal edge count)")
    axB.set_xticks(list(x))
    axB.set_xticklabels([n.replace("_", " ") for n in names], fontsize=8, rotation=15)
    axB.set_ylabel("in-sample fit (Pearson)")
    axB.set_title("B  Reasoning versus random search, in-sample fit", fontsize=11)
    axB.axhline(0, color="#333333", lw=0.6)
    axB.legend(fontsize=9)

    fig.suptitle("MMC engineer-behavior, scaled: no held-out advantage across a powered corpus",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.005,
             f"Across {d['summary']['n_proposals']} proposals and {len(ne)} novel-edge instances "
             f"({len(distinct)} distinct) over 9 runs and two conditions, none is in a model that "
             f"beats a linear baseline held-out (Panel A). The reasoning step produces "
             f"better-fitting structure than random search on the informative modules, and not on "
             f"the redundant one (Panel B). Scope: this atlas, CD4+ T, these modules.",
             ha="center", fontsize=8, style="italic", wrap=True)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
