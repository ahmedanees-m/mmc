"""Work Stream A5: render the engineer-behavior characterization figure.

Reads engineer_behavior.json and draws two panels: the held-out DE-overlap of the model
against the linear baseline per module (the model does not predict), and the proposal
characterization (coherently-argued hypotheses, split by novelty, none validated as a
held-out predictive necessity). The headline is the catch rate: the fraction of coherent
proposals the gate rejects.
"""
from __future__ import annotations

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RED, GREEN, GREY = "#C62828", "#2E7D32", "#9E9E9E"


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else "engineer_behavior.json"
    out = sys.argv[2] if len(sys.argv) > 2 else "engineer_behavior.png"
    with open(path) as f:
        d = json.load(f)

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

    # Panel B: proposal characterization
    n_prop = d["n_coherent_proposals"]
    n_novel = d["n_novel_hypotheses"]
    n_val = d["n_novel_validated_held_out"]
    catch = d.get("catch_rate")
    axB.bar([0], [n_prop], color=GREEN, width=0.5, label="coherently-argued proposals")
    axB.bar([1], [n_novel], color="#F9A825", width=0.5, label="novel (non-textbook) hypotheses")
    axB.bar([2], [n_val], color=RED, width=0.5, label="validated held-out (predictive necessity)")
    axB.set_xticks([0, 1, 2])
    axB.set_xticklabels(["proposed\n(all coherent)", "novel", "validated\nheld-out"], fontsize=9)
    axB.set_ylabel("number of structural hypotheses")
    axB.set_title("B  Plausibility does not track predictive validity", fontsize=11)
    for i, v in enumerate([n_prop, n_novel, n_val]):
        axB.text(i, v + 0.3, str(v), ha="center", fontsize=11, fontweight="bold")
    catch_txt = "n/a" if catch is None else f"{catch*100:.0f}%"
    axB.text(0.5, 0.9, f"catch rate {catch_txt}\n(coherent proposals the gate rejects)",
             transform=axB.transAxes, ha="center", fontsize=10,
             bbox=dict(boxstyle="round", fc="#FFF3E0", ec=RED))

    fig.suptitle("MMC engineer-behavior: the held-out gate supplies the calibration the "
                 "rationale lacks", fontsize=12, fontweight="bold")
    fig.text(0.5, 0.005, "Across the loop's coherently-argued structural hypotheses, none is "
             "validated as a held-out predictive necessity; the gate rejects every "
             "non-predictive proposal. Scope: this atlas, CD4+ T, these modules.",
             ha="center", fontsize=9, style="italic")
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
