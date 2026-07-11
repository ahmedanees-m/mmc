"""Streamlit presentation of MMC results.

Runs on a clean clone with Streamlit alone; it renders precomputed results and requires no
API key, data store, or GPU. Three views: interrogate the circuit, held-out evaluation, and
the limit map.

    streamlit run demo/app.py
"""
from __future__ import annotations

import json
import os

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
LIMIT_MAP = os.path.join(HERE, "..", "paper", "mmc_limit_map.png")
ENGINEER_FIG = os.path.join(HERE, "..", "paper", "engineer_behavior.png")
ENGINEER_SCALED_FIG = os.path.join(HERE, "..", "paper", "engineer_behavior_scaled.png")
DEMO_INTERV = os.path.join(HERE, "demo_interventions.json")

st.set_page_config(page_title="MMC: mechanistic discovery loop", layout="wide")

GREEN, RED = "#2E7D32", "#C62828"


def circuit_dot(edges, highlight=None):
    highlight = highlight or set()
    lines = ["digraph {", "rankdir=LR;", "bgcolor=transparent;",
             'node [shape=box style="rounded,filled" fontname="Helvetica" '
             'fillcolor="#F5F5F5" color="#BDBDBD"];',
             'edge [fontname="Helvetica" penwidth=2];']
    nodes = {g for e in edges for g in (e[0], e[1])}
    for n in sorted(nodes):
        fc = '"#FFF3E0"' if n in highlight else '"#F5F5F5"'
        col = f'"{RED}"' if n in highlight else '"#BDBDBD"'
        lines.append(f'"{n}" [fillcolor={fc} color={col}];')
    for r, t, s in edges:
        color = GREEN if s > 0 else RED
        arrow = "normal" if s > 0 else "tee"
        lines.append(f'"{r}" -> "{t}" [color="{color}" arrowhead={arrow}];')
    lines.append("}")
    return "\n".join(lines)


def intervention_dot(genes, edges, response, path_edges):
    path_set = {(r, t) for r, t in path_edges}
    lines = ["digraph {", "rankdir=LR;", "bgcolor=transparent;",
             'node [shape=box style="rounded,filled" fontname="Helvetica"];',
             'edge [fontname="Helvetica"];']
    for g in genes:
        d = response.get(g, 0.0)
        if d <= -0.5:
            fill, col = "#FFCDD2", RED
        elif d >= 0.5:
            fill, col = "#C8E6C9", GREEN
        else:
            fill, col = "#F5F5F5", "#BDBDBD"
        label = f"{g}\\n{d:+.1f}" if abs(d) >= 0.05 else g
        lines.append(f'"{g}" [label="{label}" fillcolor="{fill}" color="{col}"];')
    for r, t, s in edges:
        on = (r, t) in path_set
        color = (GREEN if s > 0 else RED) if on else "#D0D0D0"
        pen = 3 if on else 1
        arrow = "normal" if s > 0 else "tee"
        lines.append(f'"{r}" -> "{t}" [color="{color}" penwidth={pen} arrowhead={arrow}];')
    lines.append("}")
    return "\n".join(lines)


# ------------------------------ header ------------------------------
st.title("MMC: mechanistic discovery loop for T-cell gene regulation")
st.markdown(
    "MMC reads a genome-scale immune-perturbation atlas, proposes interpretable, runnable "
    "models of gene-regulatory circuits, fits them to the measured perturbation responses, and "
    "evaluates each against simple baselines on held-out data. Structure and logic are set by "
    "the reasoning step; magnitudes are set by the optimizer."
)
st.caption("Scope: the Zhu 2025 genome-scale CD4+ T-cell Perturb-seq atlas, these modules, "
           "CD4+ T cells. No prediction win or disease discovery is claimed.")
st.divider()

view = st.radio("View", ["Interrogate the circuit", "Held-out evaluation", "Limit map"],
                horizontal=True, label_visibility="collapsed")

# ------------------------------ Interrogate the circuit ------------------------------
if view == "Interrogate the circuit":
    st.header("Interrogate the circuit: the Th2 / GATA3 axis")
    st.markdown(
        "The loop reconstructs a runnable, interpretable model of the Th2 / GATA3 axis "
        "(associated with allergy, asthma, and atopic disease) from the atlas. Select an "
        "intervention below: the structural model simulates the response and traces the causal "
        "path through the validated circuit. A model that can be interrogated, not a black box."
    )
    if not os.path.exists(DEMO_INTERV):
        st.warning("demo_interventions.json not found; run scripts/demo_precompute.py")
    else:
        with open(DEMO_INTERV) as fh:
            data = json.load(fh)
        choice = st.selectbox("Intervention",
                              ["(unperturbed circuit)"] + list(data["interventions"].keys()))
        c1, c2 = st.columns([3, 2])
        if choice == "(unperturbed circuit)":
            with c1:
                st.graphviz_chart(
                    intervention_dot(data["genes"], data["edges"],
                                     {g: 0.0 for g in data["genes"]}, []),
                    use_container_width=True)
                st.caption("The validated Th2 / GATA3 circuit; every edge was confirmed "
                           "held-out-necessary by the edge-ablation control. Select an "
                           "intervention to simulate.")
            with c2:
                st.markdown("**How to read it.** Green arrow: activation. Red bar: repression. "
                            "After an intervention, each node's color and label show the "
                            "simulated log2 fold-change, and the highlighted path is the causal "
                            "route the effect travels.")
        else:
            iv = data["interventions"][choice]
            with c1:
                st.graphviz_chart(
                    intervention_dot(data["genes"], data["edges"], iv["response"], iv["path"]),
                    use_container_width=True)
                st.caption("Node label: simulated log2 fold-change. Highlighted path: the causal "
                           "route through validated edges. Red node: down, green node: up.")
            with c2:
                st.markdown(f"**Cytokine readouts, simulated ({iv['kind']}):**")
                for m, g in zip(st.columns(len(data["readouts"])), data["readouts"]):
                    m.metric(g, f"{iv['response'][g]:+.2f}")
                if iv["composed"] and iv["additive_baseline"]:
                    add = iv["additive_baseline"]
                    st.markdown("**Structural model versus the additive baseline:**")
                    st.table({"readout": data["readouts"],
                              "structural": [f"{iv['response'][g]:+.2f}"
                                             for g in data["readouts"]],
                              "additive (sum)": [f"{add[g]:+.2f}" for g in data["readouts"]]})
                    st.info("STAT6 acts through GATA3. Once GATA3 is knocked down, also knocking "
                            "down STAT6 adds little, and the structural model captures this "
                            "epistasis; the additive baseline sums the two knockdowns and "
                            "over-predicts the drop. Composing interventions is what a mechanistic "
                            "model can do and an additive baseline cannot.")
                else:
                    st.caption("The response propagates only through validated edges. A "
                               "single-gene baseline gives these numbers but cannot trace the "
                               "path or compose a second intervention.")

# ------------------------------ Held-out evaluation ------------------------------
elif view == "Held-out evaluation":
    st.header("Held-out evaluation: cytokine-production module")
    st.markdown(
        "On the cytokine-production module, the loop proposed a novel hypothesis and reasoned "
        "about it from the data. The held-out evaluation (leave-one-perturbation-out) then "
        "measures whether the model predicts better than a simple baseline on unseen data."
    )
    left, right = st.columns(2)
    with left:
        st.subheader("Proposed hypothesis")
        st.markdown("**Proposal:**")
        st.info("Use STK11 (LKB1) as a metabolic hub that gates cytokines and controls the "
                "mitochondrial/metabolic genes, which in turn each feed one cytokine.")
        st.markdown("**Revised from the residuals:**")
        st.info("The STK11 residuals consistently show that knockdown increases chemokine "
                "output, indicating STK11 as a repressor of the chemokines CCL3 / CCL4 / CXCL8.")
        st.graphviz_chart(circuit_dot(
            [("STK11", "IL2", -1), ("STK11", "IFNG", -1), ("STK11", "CXCL8", -1),
             ("STK11", "CCL3", -1), ("STK11", "CCL4", -1)], highlight={"STK11"}),
            use_container_width=True)
        st.caption("The STK11 edges are individually supported: an edge-ablation control flags "
                   "them as predictively necessary, comparable to textbook edges.")
    with right:
        st.subheader("Held-out result")
        st.markdown("Leave-one-perturbation-out, model versus simple baselines:")
        st.table({
            "method": ["**MMC model**", "linear baseline", "mean baseline"],
            "held-out DE-overlap [95% CI]": ["**0.18 [0.10, 0.27]**",
                                             "0.45 [0.33, 0.57]", "0.37 [0.26, 0.49]"],
        })
        st.error("Proposed, not certified. The STK11 edges are individually supported, but the "
                 "model built from them predicts below a linear baseline on held-out data, with "
                 "separated confidence intervals. Edge-level support does not imply predictive "
                 "advantage over the baseline.")
        st.info("The distinction the evaluation enforces: a supported marginal effect is not the "
                "same as a model that predicts better than a simple baseline.")

    st.divider()
    st.subheader("Measured across the modules")
    st.markdown(
        "The loop's novel hypotheses are individually supported: the same edge-ablation gate "
        "flags the novel STK11 edges as predictively necessary, as it flags textbook edges. "
        "Across a powered corpus (25 proposals and 76 distinct novel hypotheses over 9 runs and "
        "two conditions), none is in a model that beats a linear baseline held-out."
    )
    m1, m2, m3 = st.columns(3)
    m1.metric("Novel hypotheses validated", "0 of 76", help="Wilson 95% CI [0, 4.8%]")
    m2.metric("Module-conditions beating linear", "0 of 9", help="Wilson 95% CI [0, 30%]")
    m3.metric("Reasoning vs random (in-sample)", "0.20 vs 0.07", help="mean in-sample Pearson")
    if os.path.exists(ENGINEER_SCALED_FIG):
        st.image(ENGINEER_SCALED_FIG, use_container_width=True)
    elif os.path.exists(ENGINEER_FIG):
        st.image(ENGINEER_FIG, use_container_width=True)
    st.caption("Edge-level support is not predictive advantage over a baseline. The reasoning "
               "step produces better-fitting structure than random search on the informative "
               "modules; the held-out evaluation separates a supported marginal effect from a "
               "predictive model.")

# ------------------------------ Limit map ------------------------------
else:
    st.header("Limit map")
    st.markdown(
        "MMC maps the regimes where a mechanistic or AI model does and does not beat simple "
        "baselines on held-out data. This addresses why such models frequently fail to beat "
        "baselines: they are commonly evaluated where no method has an advantage, and in-sample "
        "fit does not imply held-out prediction."
    )
    if os.path.exists(LIMIT_MAP):
        st.image(LIMIT_MAP, use_container_width=True)
    else:
        st.warning("Limit-map figure not found at paper/mmc_limit_map.png")
    st.markdown(
        "On single-knockdown steady-state data, the mechanistic model shows no held-out "
        "advantage over simple baselines in any measured regime, whether or not it fits "
        "in-sample. The remaining source of advantage is non-additivity (combinatorial "
        "perturbations), which single-knockdown data cannot exercise; the Norman combinatorial "
        "test shows that a model fit on single perturbations does not recover it either."
    )

st.divider()
st.markdown(
    "MMC produces interpretable, testable mechanistic hypotheses from a genome-scale atlas and "
    "reports the regime within which its outputs are reliable. It supports the selection of "
    "regulatory hypotheses for follow-up and declares the boundary beyond which its predictions "
    "are not supported."
)
