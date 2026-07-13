"""Streamlit interface for MMC results.

Runs on a clean clone with Streamlit alone. It renders precomputed results and requires no API
key, data store, or GPU. Four tabs: a deterministic replay of the discovery loop, a proposed
hypothesis evaluated against a baseline, an interactive reconstruction of the Th2 circuit, and
the supporting evidence.

    streamlit run demo/app.py
"""
from __future__ import annotations

import html
import json
import math
import os

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
LIMIT_MAP = os.path.join(HERE, "..", "paper", "mmc_limit_map.png")
DEMO_INTERV = os.path.join(HERE, "demo_interventions.json")
DEMO_LOOP = os.path.join(HERE, "loop_replay.json")
LOOP_FALLBACK = os.path.join(HERE, "..", "results", "trace", "loop_replay.json")

# palette
ON, PATHC, DIM = "#2c3e50", "#e67e22", "#c9d2db"
DOWN_F, DOWN_S = "#f8d7d5", "#c0392b"
UP_F, UP_S = "#e3f1e4", "#2e7d32"
NEU_F, NEU_S = "#eef2f7", "#516072"
RANDOM_C = "#8a97a4"
LOOP_NODE = "#1a3d5c"

# fixed layout for the reconstructed circuit
POS = {
    "STAT4": (95, 78),
    "STAT6": (95, 362),
    "TBX21": (300, 78),
    "GATA3": (430, 220),
    "IL4": (785, 78),
    "IL5": (785, 220),
    "IL13": (785, 362),
}
RX, RY = 62, 32

WHY = {
    "None (baseline)":
        "Resting state. GATA3 is active and maintains IL4, IL5, and IL13.",
    "GATA3 knockdown":
        "GATA3 activates IL4, IL5, and IL13 directly. Removing it reduces all three, with the "
        "largest drop in IL5 and IL13.",
    "STAT6 knockdown":
        "STAT6 acts upstream of GATA3. Removing STAT6 alone has a small effect on the cytokines "
        "because GATA3 remains active.",
    "GATA3 + STAT6 knockdown":
        "STAT6 acts through GATA3. With GATA3 already knocked down, removing STAT6 as well adds "
        "nothing, so IL5 stays at -5.0. The additive baseline sums the two single knockdowns and "
        "predicts -5.4.",
    "GATA3 activation":
        "Activating GATA3 raises IL4 most strongly. IL5 and IL13 move less in this module.",
}


def _trim(p1, p2, off1, off2):
    """Shorten a segment so it starts and ends outside the node ellipses."""
    x1, y1 = p1
    x2, y2 = p2
    dx, dy = x2 - x1, y2 - y1
    d = math.hypot(dx, dy) or 1.0
    ux, uy = dx / d, dy / d
    return (x1 + ux * off1, y1 + uy * off1), (x2 - ux * off2, y2 - uy * off2)


def circuit_svg(genes, edges, response, path, targets):
    path_set = {(a, b) for a, b in path}
    marks = []
    for name, col in (("on", ON), ("pathc", PATHC), ("dim", DIM)):
        marks.append(
            f'<marker id="ar-{name}" markerWidth="10" markerHeight="10" refX="8" refY="4.5" '
            f'orient="auto"><path d="M0,0 L9,4.5 L0,9 z" fill="{col}"/></marker>')
        marks.append(
            f'<marker id="te-{name}" markerWidth="6" markerHeight="14" refX="1.5" refY="7" '
            f'orient="auto"><rect x="0" y="0" width="2.6" height="14" fill="{col}"/></marker>')

    body = ['<rect x="0" y="0" width="880" height="440" rx="12" fill="#f7f9fb"/>']
    for a, b, s in edges:
        if a not in POS or b not in POS:
            continue
        (sx, sy), (ex, ey) = _trim(POS[a], POS[b], RX + 4, RX + 10)
        on_path = (a, b) in path_set
        if a in targets and not on_path:
            state = "dim"
        elif on_path:
            state = "pathc"
        else:
            state = "on"
        color = {"on": ON, "pathc": PATHC, "dim": DIM}[state]
        width = 5 if on_path else (2 if state == "dim" else 3)
        shape = "ar" if s > 0 else "te"
        body.append(
            f'<line x1="{sx:.0f}" y1="{sy:.0f}" x2="{ex:.0f}" y2="{ey:.0f}" stroke="{color}" '
            f'stroke-width="{width}" marker-end="url(#{shape}-{state})"/>')

    for g in genes:
        if g not in POS:
            continue
        x, y = POS[g]
        d = response.get(g, 0.0)
        if d <= -0.5:
            fill, stroke = DOWN_F, DOWN_S
        elif d >= 0.5:
            fill, stroke = UP_F, UP_S
        else:
            fill, stroke = NEU_F, NEU_S
        sw = 5 if g in targets else 3
        body.append(
            f'<ellipse cx="{x}" cy="{y}" rx="{RX}" ry="{RY}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{sw}"/>')
        if abs(d) >= 0.05:
            body.append(
                f'<text x="{x}" y="{y - 4}" text-anchor="middle" font-size="21" font-weight="700" '
                f'fill="{stroke}" font-family="DejaVu Sans, sans-serif">{g}</text>')
            body.append(
                f'<text x="{x}" y="{y + 18}" text-anchor="middle" font-size="17" fill="{stroke}" '
                f'font-family="DejaVu Sans, sans-serif">{d:+.1f}</text>')
        else:
            body.append(
                f'<text x="{x}" y="{y + 7}" text-anchor="middle" font-size="21" font-weight="700" '
                f'fill="{stroke}" font-family="DejaVu Sans, sans-serif">{g}</text>')

    return ('<svg viewBox="0 0 880 440" xmlns="http://www.w3.org/2000/svg" width="100%">'
            f'<defs>{"".join(marks)}</defs>{"".join(body)}</svg>')


def loop_circuit_svg(genes, edges, changed):
    """A left-to-right cascade. Activating edges are arrows, repressive edges end in a bar,
    edges whose sign changed from the previous iteration are highlighted. Non-adjacent edges
    are drawn as arcs so they do not pass through intervening nodes."""
    n = len(genes)
    width, height, y = 900, 300, 168
    margin, erx, ery = 90, 52, 30
    step = (width - 2 * margin) / max(n - 1, 1)
    xs = {g: margin + k * step for k, g in enumerate(genes)}
    idx = {g: k for k, g in enumerate(genes)}

    marks = []
    for name, col in (("act", LOOP_NODE), ("chg", PATHC)):
        marks.append(
            f'<marker id="lp-{name}" markerWidth="10" markerHeight="10" refX="9" refY="4.5" '
            f'orient="auto"><path d="M0,0 L9,4.5 L0,9 z" fill="{col}"/></marker>')

    body = [f'<rect x="0" y="0" width="{width}" height="{height}" rx="12" fill="#f7f9fb"/>']
    for e in edges:
        s, d, sign = e["src"], e["dst"], e["sign"]
        if s not in xs or d not in xs:
            continue
        is_chg = (s, d) in changed
        col = PATHC if is_chg else (DOWN_S if sign < 0 else LOOP_NODE)
        w = 6 if is_chg else 3.5
        dist = idx[d] - idx[s]
        forward = xs[d] > xs[s]
        sx = xs[s] + (erx if forward else -erx)
        ex = xs[d] - (erx if forward else -erx)
        if abs(dist) == 1:
            if sign < 0:
                body.append(f'<line x1="{sx:.0f}" y1="{y}" x2="{ex:.0f}" y2="{y}" '
                            f'stroke="{col}" stroke-width="{w}"/>')
                body.append(f'<line x1="{ex:.0f}" y1="{y - 13}" x2="{ex:.0f}" y2="{y + 13}" '
                            f'stroke="{col}" stroke-width="{w}"/>')
            else:
                body.append(f'<line x1="{sx:.0f}" y1="{y}" x2="{ex:.0f}" y2="{y}" stroke="{col}" '
                            f'stroke-width="{w}" marker-end="url(#lp-{"chg" if is_chg else "act"})"/>')
            lx, ly = (sx + ex) / 2, y - 20
        else:
            up = dist > 0
            cx = (xs[s] + xs[d]) / 2
            cy = y - (46 + 16 * abs(dist)) if up else y + (46 + 16 * abs(dist))
            sy = y - ery * 0.5 if up else y + ery * 0.5
            sxa = xs[s] + (erx * 0.5 if forward else -erx * 0.5)
            exa = xs[d] - (erx * 0.5 if forward else -erx * 0.5)
            path = f'M {sxa:.0f},{sy:.0f} Q {cx:.0f},{cy:.0f} {exa:.0f},{sy:.0f}'
            if sign < 0:
                body.append(f'<path d="{path}" fill="none" stroke="{col}" stroke-width="{w}"/>')
                body.append(f'<line x1="{exa:.0f}" y1="{sy - 11:.0f}" x2="{exa:.0f}" '
                            f'y2="{sy + 11:.0f}" stroke="{col}" stroke-width="{w}"/>')
            else:
                body.append(f'<path d="{path}" fill="none" stroke="{col}" stroke-width="{w}" '
                            f'marker-end="url(#lp-{"chg" if is_chg else "act"})"/>')
            lx, ly = cx, (cy - 8 if up else cy + 18)
        if is_chg:
            body.append(f'<text x="{lx:.0f}" y="{ly:.0f}" text-anchor="middle" font-size="16" '
                        f'font-weight="700" fill="{PATHC}" font-family="DejaVu Sans, sans-serif">'
                        f'flipped</text>')

    for g in genes:
        x = xs[g]
        body.append(f'<ellipse cx="{x:.0f}" cy="{y}" rx="{erx}" ry="{ery}" fill="#eaf2fb" '
                    f'stroke="{LOOP_NODE}" stroke-width="3.5"/>')
        body.append(f'<text x="{x:.0f}" y="{y + 6}" text-anchor="middle" font-size="17" '
                    f'font-weight="700" fill="{LOOP_NODE}" font-family="DejaVu Sans, sans-serif">'
                    f'{g}</text>')
    return (f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
            f'width="100%"><defs>{"".join(marks)}</defs>{"".join(body)}</svg>')


_RES_LABEL = {"wrong_sign": "sign is wrong", "missing_effect": "effect missing",
              "spurious_effect": "spurious effect"}


def residuals_html(residuals):
    rows = []
    for r in residuals:
        pat = r.get("pattern", "")
        cls = "res-mag" if pat == "magnitude" else "res-sign"
        tag = _RES_LABEL.get(pat, pat.replace("_", " "))
        rows.append(
            f'<div class="{cls}">{r["perturbation"]} knockdown, {r["target"]}: '
            f'model {r["predicted"]:+.2f}, observed {r["observed"]:+.2f} ({tag})</div>')
    return "".join(rows)


def bars_svg(rows, vmax, width=460, rowh=52, gutter=170):
    height = rowh * len(rows) + 10
    body = [f'<rect x="0" y="0" width="{width}" height="{height}" rx="8" fill="#f7f9fb"/>']
    for i, (label, val, color) in enumerate(rows):
        yy = 8 + i * rowh
        bw = max(2.0, (val / vmax) * (width - gutter - 60))
        body.append(
            f'<text x="0" y="{yy + 24}" font-size="17" fill="#3a4652" '
            f'font-family="DejaVu Sans, sans-serif">{label}</text>')
        body.append(
            f'<rect x="{gutter}" y="{yy + 8}" width="{bw:.0f}" height="26" rx="5" fill="{color}"/>')
        body.append(
            f'<text x="{gutter + bw + 8:.0f}" y="{yy + 27}" font-size="17" font-weight="700" '
            f'fill="{color}" font-family="DejaVu Sans, sans-serif">{val:.2f}</text>')
    return (f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" width="100%">'
            f'{"".join(body)}</svg>')


def load_json(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return None


st.set_page_config(page_title="MMC", layout="wide")
st.markdown(
    """
<style>
  #MainMenu, footer, header {visibility: hidden;}
  .block-container {padding-top: 2rem; max-width: 1180px;}
  html, body, [class*="css"] {font-size: 18px;}
  [data-testid="stMetricValue"] {font-size: 2.5rem; font-weight: 700;}
  [data-testid="stMetricLabel"] {font-size: 1rem; color: #4a5568;}
  .mmc-title {font-size: 2.1rem; font-weight: 800; margin-bottom: 0.1rem;}
  .mmc-sub {font-size: 1.08rem; color: #5a6572; margin-top: 0; margin-bottom: 0.4rem;}
  .iter {font-size: 1.35rem; font-weight: 800; margin: 0.3rem 0 0.2rem;}
  .note {font-size: 1.05rem; color: #2b3640; background: #f4f6f8; padding: 14px 16px;
         border-radius: 8px; border-left: 5px solid #37474f; margin-top: 10px;}
  .think {font-size: 1.1rem; line-height: 1.5; color: #43340f; background: #fffdf5;
          padding: 16px 18px; border-radius: 8px; border-left: 5px solid #e67e22;
          min-height: 120px;}
  .edit {font-size: 1.05rem; font-weight: 600; color: #14385c; background: #eef6ff;
         padding: 12px 16px; border-radius: 8px; border-left: 5px solid #2c6fb5;
         margin-top: 10px;}
  .res-sign {font-family: DejaVu Sans Mono, monospace; font-size: 1rem; color: #b21f12;
             font-weight: 600; line-height: 1.7;}
  .res-mag {font-family: DejaVu Sans Mono, monospace; font-size: 1rem; color: #8a8f94;
            line-height: 1.7;}
  .flag-red {background: #fdecea; color: #8e1b12; border-left: 5px solid #c0392b;
             padding: 12px 16px; border-radius: 8px; font-weight: 600; margin-top: 10px;}
  .flag-green {background: #e8f5e9; color: #1b5e20; border-left: 5px solid #2e7d32;
               padding: 12px 16px; border-radius: 8px; font-weight: 600; margin-top: 10px;}
</style>
""",
    unsafe_allow_html=True,
)

st.markdown('<p class="mmc-title">MMC: a mechanistic model-discovery loop for T-cell regulation</p>',
            unsafe_allow_html=True)
st.markdown('<p class="mmc-sub">Reconstructed from the Zhu 2025 CD4+ T-cell Perturb-seq atlas. '
            'The reasoning step proposes and revises structure; an optimizer sets magnitudes; a '
            'held-out gate decides what is certified.</p>', unsafe_allow_html=True)

data = load_json(DEMO_INTERV)
loop = load_json(DEMO_LOOP) or load_json(LOOP_FALLBACK)

tab_loop, tab_stk, tab_circuit, tab_evidence = st.tabs(
    ["Discovery loop", "STK11 proposal", "Circuit", "Evidence"])

# ------------------------------ Discovery loop ------------------------------
with tab_loop:
    st.header("Discovery loop")
    st.caption("Deterministic replay of a captured run on the TCR signalosome core in Stim8hr. "
               "The reasoning step proposes a structure, the simulator fits it, the fit-versus-"
               "structure gate surfaces the structural residuals, and the reasoning step revises "
               "the structure. Training data only.")
    if not loop:
        st.info("No captured trace found. Run scripts/capture_loop_trace.py (with MMC_ZHU_STORE "
                "and ANTHROPIC_API_KEY) to write results/trace/loop_replay.json.")
    else:
        iters = loop["iterations"]
        if "li" not in st.session_state:
            st.session_state.li = 0
        i = min(st.session_state.li, len(iters) - 1)
        it = iters[i]

        changed = set()
        if i > 0:
            prev = {(e["src"], e["dst"]): e["sign"] for e in iters[i - 1]["edges"]}
            for e in it["edges"]:
                if prev.get((e["src"], e["dst"]), e["sign"]) != e["sign"]:
                    changed.add((e["src"], e["dst"]))

        st.markdown(f'<p class="iter">Iteration {it["n"]} of {len(iters)}</p>',
                    unsafe_allow_html=True)
        st.markdown(loop_circuit_svg(it["genes"], it["edges"], changed), unsafe_allow_html=True)

        left, right = st.columns(2, gap="large")
        with left:
            st.markdown("##### What the fit revealed (residuals shown to the reasoning step)")
            if it["residuals"]:
                st.markdown(f'<div class="note">{residuals_html(it["residuals"])}</div>',
                            unsafe_allow_html=True)
            else:
                st.markdown('<div class="flag-green">No structural residuals remain.</div>',
                            unsafe_allow_html=True)
            mm = st.columns(2)
            prev_fit = iters[i - 1]["fit"] if i > 0 else None
            prev_ns = iters[i - 1]["n_structural"] if i > 0 else None
            mm[0].metric("Training loss", f'{it["fit"]:.2f}',
                         delta=None if prev_fit is None else f'{it["fit"] - prev_fit:+.2f}',
                         delta_color="inverse")
            mm[1].metric("Structural residuals", it["n_structural"],
                         delta=None if prev_ns is None else f'{it["n_structural"] - prev_ns:+d}',
                         delta_color="inverse")
        with right:
            st.markdown("##### The reasoning step")
            st.markdown(f'<div class="think">{html.escape(it["rationale"])}</div>',
                        unsafe_allow_html=True)
            if it.get("edit"):
                e = it["edit"]
                parts = []
                if isinstance(e, dict):
                    if e.get("removed"):
                        parts.append("removed " + ", ".join(e["removed"]))
                    if e.get("added"):
                        parts.append("added " + ", ".join(e["added"]))
                    label = "; ".join(parts) if parts else "structure revised"
                else:
                    label = str(e)
                st.markdown(f'<div class="edit">Structural edit: {html.escape(label)}</div>',
                            unsafe_allow_html=True)

        st.progress((i + 1) / len(iters))
        nav = st.columns([1, 1, 4])
        if nav[0].button("Next iteration", type="primary", disabled=i >= len(iters) - 1):
            st.session_state.li = i + 1
            st.rerun()
        if nav[1].button("Restart"):
            st.session_state.li = 0
            st.rerun()
        if i >= len(iters) - 1:
            st.markdown('<div class="note">The structure converged. What the loop set is the '
                        'structure and the logic form; the optimizer set the magnitudes. This '
                        'structure is what the held-out gate evaluates next.</div>',
                        unsafe_allow_html=True)

# ------------------------------ STK11 proposal ------------------------------
with tab_stk:
    st.markdown("#### A proposed hypothesis, evaluated against a baseline")
    left, right = st.columns(2, gap="large")
    with left:
        st.markdown("**Model proposal (verbatim)**")
        st.info("Use STK11 (LKB1) as a metabolic hub that both gates cytokines and controls the "
                "mitochondrial and metabolic genes, which in turn each feed one cytokine.")
        st.info("The STK11 residuals consistently show that knockdown increases chemokine output, "
                "indicating STK11 as a repressor of the chemokines CCL3, CCL4, and CXCL8.")
        st.markdown('<div class="note"><b>Edge support.</b> An edge-ablation control flags the '
                    'STK11 edges as predictively necessary, like the textbook edges. Edge-level '
                    'support is distinct from predictive advantage.</div>', unsafe_allow_html=True)
    with right:
        st.markdown("**Held-out evaluation (DE-overlap, higher is better)**")
        st.markdown(
            bars_svg([("MMC model", 0.18, DOWN_S), ("Linear baseline", 0.45, ON)], 0.5),
            unsafe_allow_html=True)
        mm = st.columns(2)
        mm[0].metric("MMC model (held-out)", "0.18")
        mm[1].metric("Linear baseline", "0.45")
        st.markdown('<div class="flag-red">Not certified. The model predicts below a linear '
                    'baseline on held-out data, with separated confidence intervals.</div>',
                    unsafe_allow_html=True)

# ------------------------------ Circuit ------------------------------
with tab_circuit:
    st.caption("The reconstructed Th2 / GATA3 circuit. Select an intervention to simulate the "
               "response and trace the causal path. Validated edges only.")
    if not data:
        st.warning("demo_interventions.json not found.")
    else:
        avail = set(data["interventions"])
        iv_order = [("None (baseline)", None)] + [
            (k, k) for k in ("GATA3 knockdown", "STAT6 knockdown", "GATA3 + STAT6 knockdown",
                             "GATA3 activation") if k in avail
        ]
        iv_map = dict(iv_order)
        choice = st.radio("Intervention", [lbl for lbl, _ in iv_order], horizontal=True)
        key = iv_map[choice]
        if key is None:
            iv = {"response": {g: 0.0 for g in data["genes"]}, "path": [], "additive": None,
                  "targets": []}
        else:
            rec = data["interventions"][key]
            iv = {"response": rec["response"], "path": rec["path"],
                  "additive": rec.get("additive_baseline"), "targets": rec.get("targets", [])}

        left, right = st.columns([1.25, 1], gap="large")
        with left:
            st.markdown(
                circuit_svg(data["genes"], data["edges"], iv["response"], iv["path"],
                            iv["targets"]),
                unsafe_allow_html=True)
            st.caption("Node color: simulated log2 fold-change (red down, green up). Highlighted "
                       "edges: the causal path. Thick border marks the perturbed gene.")
        with right:
            st.markdown("#### Predicted cytokine response")
            cc = st.columns(3)
            cc[0].metric("IL5", f"{iv['response']['IL5']:+.1f}")
            cc[1].metric("IL13", f"{iv['response']['IL13']:+.1f}")
            cc[2].metric("IL4", f"{iv['response']['IL4']:+.1f}")
            add = iv["additive"]
            if add is not None:
                st.markdown("#### IL5: MMC compared with the additive baseline")
                mm = st.columns(2)
                mmc_il5 = iv["response"]["IL5"]
                add_il5 = add["IL5"]
                mm[0].metric("MMC (mechanistic)", f"{mmc_il5:+.1f}")
                diff = add_il5 - mmc_il5
                mm[1].metric("Additive baseline", f"{add_il5:+.1f}",
                             delta=f"{diff:+.1f} vs MMC" if abs(diff) >= 0.05 else "matches MMC",
                             delta_color="off")
            st.markdown(f'<div class="note"><b>Why:</b> {WHY[choice]}</div>',
                        unsafe_allow_html=True)
            if key == "GATA3 + STAT6 knockdown":
                st.markdown('<div class="flag-green">MMC holds IL5 at -5.0. The additive baseline '
                            'predicts -5.4, over-predicting the effect by 0.4.</div>',
                            unsafe_allow_html=True)

# ------------------------------ Evidence ------------------------------
with tab_evidence:
    st.markdown("#### Across a powered corpus")
    cc = st.columns(2)
    cc[0].metric("Novel hypotheses that beat the baseline (held-out)", "0 of 76")
    cc[1].metric("Module-conditions that beat the linear baseline", "0 of 9")
    st.caption("Wilson 95% confidence interval [0, 4.8%], across 25 proposals over 9 runs and two "
               "conditions.")

    left, right = st.columns(2, gap="large")
    with left:
        st.markdown("**Reasoning compared with random search (in-sample Pearson)**")
        st.markdown(
            bars_svg([("Proposed structure", 0.20, ON), ("Random, equal edges", 0.07, RANDOM_C)],
                     0.25),
            unsafe_allow_html=True)
        st.caption("The proposed structure fits in-sample better than a random structure of equal "
                   "edge count. This checks the proposal step; it is not a held-out prediction "
                   "result.")
    with right:
        if os.path.exists(LIMIT_MAP):
            st.image(LIMIT_MAP, width="stretch")
            st.caption("Limit map: where mechanistic and AI models beat simple baselines on "
                       "held-out data and where they do not.")
        else:
            st.warning("Limit-map figure not found at paper/mmc_limit_map.png")

st.divider()
st.caption("Precomputed from the structural backend. Interventions restricted to ablation-validated "
           "edges. Scope: the Zhu 2025 CD4+ T-cell Perturb-seq atlas, Th2 module and TCR core.")
