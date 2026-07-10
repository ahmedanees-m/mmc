"""MMC demo surface: the honest AI biological engineer, in three beats.

Runs on a clean clone with only streamlit (it presents captured results; no API key,
store, or GPU needed). The three beats: Claude builds a runnable disease circuit, Claude
catches its own plausible-but-wrong hypothesis, and the tool declares its own limits.

    streamlit run demo/app.py

Framing rule (non-negotiable): no claimed prediction win, no claimed disease discovery.
The hero is autonomous interpretable modeling plus self-correction plus honest limits.
"""
from __future__ import annotations

import os

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
LIMIT_MAP = os.path.join(HERE, "..", "paper", "mmc_limit_map.png")
ENGINEER_FIG = os.path.join(HERE, "..", "paper", "engineer_behavior.png")

st.set_page_config(page_title="MMC, the honest AI biological engineer", layout="wide")

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
st.title("MMC: an AI biological engineer that knows when it is wrong")
st.markdown(
    "It reads the newest immune-disease atlas, **proposes interpretable, runnable models** "
    "of disease-driving regulatory circuits, **tests them against the real data**, and "
    "**refuses to certify the hypotheses it cannot stand behind**, catching its own "
    "plausible-but-wrong ideas. The first mechanistic-hypothesis engine a target team can "
    "actually trust, *because it tells you when not to.*"
)
st.caption("Scope: the Zhu 2025 genome-scale CD4+ T-cell Perturb-seq atlas; these modules; "
           "CD4+ T cells. No prediction win or disease discovery is claimed.")
st.divider()

beat = st.radio("Demo", ["1 · It builds", "2 · It holds itself to the baseline",
                         "3 · It knows its limits"],
                horizontal=True, label_visibility="collapsed")

# ------------------------------ Beat 1 ------------------------------
if beat.startswith("1"):
    st.header("Beat 1 · Watch Claude do science")
    st.markdown(
        "Given only the atlas, Claude autonomously writes a **runnable, interpretable "
        "model** of the **Th2 / GATA3 axis**, the circuit driving allergy, asthma, and "
        "atopic disease, as signed edges and logic it can simulate, then fits it, reads "
        "the residuals, and revises."
    )
    th2_edges = [
        ("STAT6", "GATA3", 1), ("TBX21", "GATA3", -1), ("GATA3", "TBX21", -1),
        ("STAT4", "TBX21", 1), ("TBX21", "STAT4", 1), ("IL4", "STAT6", 1),
        ("GATA3", "IL4", 1), ("GATA3", "IL5", 1), ("GATA3", "IL13", 1),
    ]
    c1, c2 = st.columns([3, 2])
    with c1:
        st.graphviz_chart(circuit_dot(th2_edges), use_container_width=True)
        st.caption("Green arrow = activation, red bar = repression. Claude proposed this "
                   "structure; the optimizer set the magnitudes.")
    with c2:
        st.markdown("**Claude's rationale (verbatim):**")
        st.info("The circuit is centered on GATA3 as the Th2 master regulator with mutual "
                "antagonism against the TBX21/STAT4 Th1 axis and a feed-forward "
                "IL4 to STAT6 to GATA3 loop, using product-gate terms only where "
                "activator/repressor logic is non-monotone.")
        st.metric("In-sample fit of the circuit", "Pearson 0.93",
                  help="How well the fitted circuit reproduces the training knockdowns. "
                       "This is a fit, not a prediction claim.")
        st.caption("The GATA3 to IL4/IL5/IL13 amplification and the GATA3 to TBX21 mutual "
                   "antagonism are textbook Th2 biology, recovered from the data.")

# ------------------------------ Beat 2 ------------------------------
elif beat.startswith("2"):
    st.header("Beat 2 · It holds itself to the baseline")
    st.markdown(
        "On the cytokine-production module, Claude proposed a **novel, disease-relevant** "
        "hypothesis and reasoned about it from the data. Then the held-out gate was asked the "
        "only question that matters: **does the model beat a simple baseline on data it never "
        "saw?**"
    )
    left, right = st.columns(2)
    with left:
        st.subheader("Claude proposes, and reasons from the data")
        st.markdown("**Proposal:**")
        st.info("...use **STK11 (LKB1)** as a metabolic hub that both gates cytokines and "
                "controls the mitochondrial/metabolic genes, which in turn each feed one "
                "cytokine.")
        st.markdown("**Then it reads the residuals and revises:**")
        st.info("The STK11 residuals consistently show that knockdown *increases* chemokine "
                "output, so **STK11 acts as a repressor of the chemokines CCL3 / CCL4 / "
                "CXCL8**.")
        st.graphviz_chart(circuit_dot(
            [("STK11", "IL2", -1), ("STK11", "IFNG", -1), ("STK11", "CXCL8", -1),
             ("STK11", "CCL3", -1), ("STK11", "CCL4", -1)], highlight={"STK11"}),
            use_container_width=True)
        st.caption("A coherent, disease-relevant hypothesis, and the STK11 edges are **real**: "
                   "an edge-ablation control flags them predictively necessary, just like "
                   "textbook edges. Not a hallucination.")
    with right:
        st.subheader("The gate refuses to certify the model")
        st.markdown("Held-out prediction (leave-one-perturbation-out): **does the mechanistic "
                    "model beat a simple baseline on data it never saw?**")
        st.table({
            "method": ["**MMC model**", "linear baseline", "mean baseline"],
            "held-out DE-overlap [95% CI]": ["**0.18 [0.10, 0.27]**",
                                             "0.45 [0.33, 0.57]", "0.37 [0.26, 0.49]"],
        })
        st.error("**PROPOSED, NOT CERTIFIED.** The STK11 edges are individually grounded, but "
                 "the mechanistic model built from them predicts *worse* than a simple linear "
                 "baseline, with cleanly separated confidence intervals. Grounded mechanism is "
                 "not predictive advantage, so the engine will not send a target team chasing "
                 "it as a discovery.")
        st.success("An AI scientist that will not overclaim even its own grounded hypotheses, "
                   "the exact thing the field is afraid AI cannot do. The honesty is the "
                   "result.")

    st.divider()
    st.subheader("This was not a one-off: we measured it")
    st.markdown(
        "Across the modules, the loop's novel hypotheses are **individually grounded**, not "
        "hallucinated: the same edge-ablation gate flags the novel STK11 edges predictively "
        "necessary, exactly as it flags textbook edges. Yet the mechanistic model still does "
        "**not beat a linear baseline** held-out on any module, and the module-level gate is "
        "the calibration that reveals the gap."
    )
    m1, m2, m3 = st.columns(3)
    m1.metric("STK11 edges flagged necessary", "3 of 3")
    m2.metric("Modules beating linear held-out", "0 of 2")
    m3.metric("Cytokine model vs linear", "0.18 vs 0.45")
    if os.path.exists(ENGINEER_FIG):
        st.image(ENGINEER_FIG, use_container_width=True)
    st.caption("Plausible, edge-grounded mechanism is not predictive advantage over a "
               "baseline. The held-out gate is what separates a real marginal effect from a "
               "model worth trusting.")

# ------------------------------ Beat 3 ------------------------------
else:
    st.header("Beat 3 · It knows its limits: trust as rigor")
    st.markdown(
        "MMC does not just fail quietly; it **maps where mechanism can and cannot be "
        "trusted**, and says so. This resolves a live confusion in the field (why do "
        "mechanistic and foundation models keep failing to beat simple baselines?): "
        "they are usually tested where nothing can win, **and** fitting is not predicting."
    )
    if os.path.exists(LIMIT_MAP):
        st.image(LIMIT_MAP, use_container_width=True)
    else:
        st.warning("limit-map figure not found at paper/mmc_limit_map.png")
    st.markdown(
        "**The finding:** on single-knockdown steady-state data, mechanism has no held-out "
        "advantage over simple baselines **in any regime**, whether or not it can fit. "
        "The only headroom is genuine non-additivity (combinatorial perturbations), which "
        "single-knockdown data structurally cannot reach. A tool you can trust *because* it "
        "declares its own boundary."
    )

st.divider()
st.markdown(
    "**For a target team:** interpretable, testable mechanistic hypotheses from the newest "
    "atlas, and an AI that will not send you chasing plausible-but-wrong targets. Most "
    "drugs fail because the target was wrong; MMC attacks that, on data you can trust, with "
    "a model that tells you when not to. **Trustworthy mechanistic discovery, by construction.**"
)
