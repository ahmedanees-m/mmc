"""Render the engineer-behavior characterization figure (corrected, two honest panels).

Panel A: held-out DE-overlap of the model against the linear baseline, per module -- the
mechanistic model does not beat a simple baseline held-out. Panel B: the edge-ablation
positive control -- the same held-out gate flags the loop's novel hypotheses (STK11 ->
chemokine) as predictively necessary, exactly as it flags textbook edges, so the novel
hypotheses are individually grounded, not hallucinated. The corrected message: plausible,
edge-grounded mechanism does not compose into a model that beats a baseline, and the
module-level held-out gate is the calibration that reveals the gap.

Reads engineer_behavior.json (Panel A) and gate_discrimination.json (Panel B).
"""
from __future__ import annotations

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RED, GREEN, GREY, ORANGE = "#C62828", "#2E7D32", "#9E9E9E", "#F9A825"


def main() -> None:
    eb = sys.argv[1] if len(sys.argv) > 1 else "engineer_behavior.json"
    gd = sys.argv[2] if len(sys.argv) > 2 else "gate_discrimination.json"
    out = sys.argv[3] if len(sys.argv) > 3 else "engineer_behavior.png"
    with open(eb) as f:
        d = json.load(f)
    with open(gd) as f:
        g = json.load(f)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: held-out DE-overlap, model vs linear, per module
    mods = d["modules"]
    labels = [m["module"].replace("_", " ") for m in mods]
    model = [m["held_out_DE_overlap"]["model"] or 0.0 for m in mods]
    linear = [m["held_out_DE_overlap"]["linear"] or 0.0 for m in mods]
    x = range(len(mods))
    axA.bar([i - 0.2 for i in x], model, width=0.4, color=RED, label="MMC model")
    axA.bar([i + 0.2 for i in x], linear, width=0.4, color=GREY, label="linear baseline")
    axA.set_xticks(list(x))
    axA.set_xticklabels(labels, fontsize=9)
    axA.set_ylabel("held-out DE-overlap (LOO)")
    axA.set_title("A  The model does not beat the linear baseline held-out", fontsize=11)
    axA.legend(fontsize=9)

    # Panel B: edge-ablation positive control (held-out ACC_DEG drop when edge removed)
    edges = sorted(g["edges"], key=lambda r: r["drop"])
    names = [e["edge"] for e in edges]
    drops = [e["drop"] for e in edges]
    colors = [GREEN if e["class"] == "textbook" else ORANGE for e in edges]
    y = range(len(edges))
    axB.barh(list(y), drops, color=colors)
    axB.axvline(g["required_drop_threshold"], color=GREY, ls="--", lw=1)
    axB.set_yticks(list(y))
    axB.set_yticklabels(names, fontsize=8)
    axB.set_xlabel("held-out ACC_DEG drop when edge removed (required if > dashed line)")
    axB.set_title("B  Novel hypotheses are individually grounded, not hallucinated",
                  fontsize=11)
    axB.plot([], [], color=GREEN, lw=6, label="textbook edge")
    axB.plot([], [], color=ORANGE, lw=6, label="novel (STK11) edge")
    axB.legend(fontsize=9, loc="lower right")

    fig.suptitle("MMC engineer-behavior: grounded mechanism does not beat a simple baseline",
                 fontsize=12, fontweight="bold")
    caption = (
        "Panel B: the edge-ablation gate flags the novel STK11 edges predictively necessary, "
        "like the textbook edges,\nso they are real marginal effects, not hallucinations. Yet "
        "the model built from them does not beat a linear\nbaseline held-out (Panel A). The "
        "module-level held-out gate is the calibration. Scope: this atlas, CD4+ T.")
    fig.text(0.5, 0.02, caption, ha="center", fontsize=8, style="italic")
    fig.tight_layout(rect=[0, 0.12, 1, 0.96])
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
