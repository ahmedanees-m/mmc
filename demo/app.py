"""Streamlit presentation of MMC results, arranged for a short screen-capture walkthrough.

Runs on a clean clone with Streamlit alone; it renders precomputed results and requires no API
key, data store, or GPU. Three screens (sidebar): the interrogable circuit, the STK11 catch, and
the evidence.

    streamlit run demo/app.py
"""
from __future__ import annotations

import json
import os

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
LIMIT_MAP = os.path.join(HERE, "..", "paper", "mmc_limit_map.png")
SCALED_FIG = os.path.join(HERE, "..", "paper", "engineer_behavior_scaled.png")
DEMO_INTERV = os.path.join(HERE, "demo_interventions.json")

GREEN, RED = "#2E7D32", "#C62828"

st.set_page_config(page_title="MMC", layout="wide")
st.markdown(
    "<style>"
    "html {font-size: 18px;}"
    "#MainMenu, footer {visibility: hidden;}"
    "[data-testid='stMetricValue'] {font-size: 2.6rem; font-weight: 700;}"
    "[data-testid='stMetricLabel'] {font-size: 1.05rem;}"
    "</style>",
    unsafe_allow_html=True,
)

# short filming labels mapped to the intervention keys in demo_interventions.json
IV_LABELS = {
    "GATA3 KD": "GATA3 knockdown",
    "STAT6 KD": "STAT6 knockdown",
    "GATA3 + STAT6 KD": "GATA3 + STAT6 knockdown",
    "GATA3 activation": "GATA3 activation",
}


def intervention_dot(genes, edges, response, path_edges):
    path_set = {(r, t) for r, t in path_edges}
    lines = ["digraph {", "rankdir=LR;", "bgcolor=transparent;", "ranksep=0.6; nodesep=0.4;",
             'node [shape=box style="rounded,filled" fontname="Helvetica" fontsize=16];',
             'edge [fontname="Helvetica" penwidth=2];']
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
        pen = 4 if on else 1
        arrow = "normal" if s > 0 else "tee"
        lines.append(f'"{r}" -> "{t}" [color="{color}" penwidth={pen} arrowhead={arrow}];')
    lines.append("}")
    return "\n".join(lines)


st.title("MMC")
st.caption("An interrogable, self-checking model of T-cell gene regulation. Scope: the Zhu 2025 "
           "CD4+ T-cell Perturb-seq atlas, these modules. No prediction win or disease discovery "
           "is claimed.")

tab1, tab2, tab3 = st.tabs(
    ["1 · Interrogate the circuit", "2 · The STK11 catch", "3 · The evidence"])

# ------------------------------ Screen 1 ------------------------------
with tab1:
    st.header("Interrogate the circuit")
    st.caption("The reconstructed Th2 / GATA3 circuit. Choose an intervention: the structural "
               "model simulates the response and traces the causal path. Validated edges only.")
    with open(DEMO_INTERV) as fh:
        data = json.load(fh)
    labels = [k for k in IV_LABELS if IV_LABELS[k] in data["interventions"]]
    choice = st.radio("Intervention", labels, horizontal=True)
    iv = data["interventions"][IV_LABELS[choice]]

    left, right = st.columns([3, 2])
    with left:
        st.graphviz_chart(
            intervention_dot(data["genes"], data["edges"], iv["response"], iv["path"]),
            use_container_width=True)
        st.caption("Node label: simulated log2 fold-change. Highlighted path: the causal route "
                   "through validated edges. Red node: down, green node: up.")
    with right:
        mmc = iv["response"]["IL5"]
        add = iv.get("additive_baseline")
        add_il5 = add["IL5"] if add else None
        st.metric("IL5  ·  MMC (mechanistic)", f"{mmc:+.1f}")
        if add_il5 is None:
            st.metric("IL5  ·  additive baseline", "n/a")
        else:
            diff = add_il5 - mmc
            st.metric("IL5  ·  additive baseline", f"{add_il5:+.1f}",
                      delta=f"{diff:+.1f}" if abs(diff) >= 0.05 else None)
        st.divider()
        cc = st.columns(2)
        cc[0].metric("IL13", f"{iv['response']['IL13']:+.1f}")
        cc[1].metric("IL4", f"{iv['response']['IL4']:+.1f}")

    if iv["composed"] and add is not None:
        st.info("STAT6 acts through GATA3. Knocking down STAT6 on top of GATA3 leaves IL5 "
                "unchanged, because GATA3 is already off; the additive baseline sums the two "
                "knockdowns and reaches -5.4. The mechanistic model composes interventions; the "
                "additive baseline cannot.")

# ------------------------------ Screen 2 ------------------------------
with tab2:
    st.header("The STK11 catch")
    left, right = st.columns(2)
    with left:
        st.subheader("Claude's proposal (verbatim)")
        st.info("Use STK11 (LKB1) as a metabolic hub that both gates cytokines and controls the "
                "mitochondrial and metabolic genes, which in turn each feed one cytokine.")
        st.info("The STK11 residuals consistently show that knockdown increases chemokine "
                "output, indicating STK11 as a repressor of the chemokines CCL3 / CCL4 / CXCL8.")
    with right:
        st.subheader("The held-out gate")
        m = st.columns(2)
        m[0].metric("MMC model, held-out", "0.18")
        m[1].metric("Linear baseline", "0.45")
        st.markdown("<h2 style='color:#C62828;margin:0.2rem 0'>NOT CERTIFIED</h2>",
                    unsafe_allow_html=True)
        st.error("The model predicts below a linear baseline on held-out data (DE-overlap 0.18 "
                 "versus 0.45), with separated confidence intervals.")
    st.markdown("**Grounded, not hallucinated.** An edge-ablation control flags the STK11 edges "
                "as predictively necessary, exactly like textbook edges. The model built from "
                "them still does not beat the baseline. Edge-level support is not predictive "
                "advantage; only the held-out comparison separates them.")

# ------------------------------ Screen 3 ------------------------------
with tab3:
    st.header("The evidence")
    left, right = st.columns([1, 1])
    with left:
        st.metric("Novel hypotheses that beat baseline held-out", "0 of 76")
        st.caption("Wilson 95% confidence interval [0, 4.8%], across 25 proposals over 9 runs "
                   "and two conditions.")
        st.metric("Module-conditions that beat linear held-out", "0 of 9")
        st.caption("Reasoning versus random search: the proposed structure fits in-sample "
                   "better than a random structure of equal edge count (mean Pearson 0.20 "
                   "versus 0.07).")
    with right:
        if os.path.exists(LIMIT_MAP):
            st.image(LIMIT_MAP, use_container_width=True)
        else:
            st.warning("Limit-map figure not found at paper/mmc_limit_map.png")
    st.caption("The limit map: where mechanistic and AI models beat simple baselines on held-out "
               "data, and where they do not. Single-perturbation prediction shows no advantage "
               "in any measured regime, and the combinatorial regime does not rescue it.")
