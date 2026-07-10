"""Streamlit presentation of MMC results.

Runs on a clean clone with Streamlit alone; it renders precomputed results and requires no
API key, data store, or GPU. Three views: circuit reconstruction, held-out evaluation, and
the limit map.

    streamlit run demo/app.py
"""
from __future__ import annotations

import os

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
LIMIT_MAP = os.path.join(HERE, "..", "paper", "mmc_limit_map.png")
ENGINEER_FIG = os.path.join(HERE, "..", "paper", "engineer_behavior.png")
ENGINEER_SCALED_FIG = os.path.join(HERE, "..", "paper", "engineer_behavior_scaled.png")

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

view = st.radio("View", ["Circuit reconstruction", "Held-out evaluation", "Limit map"],
                horizontal=True, label_visibility="collapsed")

# ------------------------------ Circuit reconstruction ------------------------------
if view == "Circuit reconstruction":
    st.header("Circuit reconstruction: the Th2 / GATA3 axis")
    st.markdown(
        "From the atlas alone, the loop writes a runnable, interpretable model of the "
        "Th2 / GATA3 axis (associated with allergy, asthma, and atopic disease) as signed edges "
        "and logic, fits it, reads the residuals, and revises the structure."
    )
    th2_edges = [
        ("STAT6", "GATA3", 1), ("TBX21", "GATA3", -1), ("GATA3", "TBX21", -1),
        ("STAT4", "TBX21", 1), ("TBX21", "STAT4", 1), ("IL4", "STAT6", 1),
        ("GATA3", "IL4", 1), ("GATA3", "IL5", 1), ("GATA3", "IL13", 1),
    ]
    c1, c2 = st.columns([3, 2])
    with c1:
        st.graphviz_chart(circuit_dot(th2_edges), use_container_width=True)
        st.caption("Green arrow: activation. Red bar: repression. The reasoning step proposed "
                   "the structure; the optimizer set the magnitudes.")
    with c2:
        st.markdown("**Proposed rationale (verbatim):**")
        st.info("The circuit is centered on GATA3 as the Th2 master regulator with mutual "
                "antagonism against the TBX21/STAT4 Th1 axis and a feed-forward IL4 to STAT6 to "
                "GATA3 loop, using product-gate terms only where activator/repressor logic is "
                "non-monotone.")
        st.metric("In-sample fit", "Pearson 0.93",
                  help="Agreement between the fitted circuit and the training knockdowns. "
                       "In-sample fit, not a prediction claim.")
        st.caption("The GATA3 to IL4/IL5/IL13 amplification and the GATA3 to TBX21 mutual "
                   "antagonism are established Th2 biology, recovered from the data.")

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
