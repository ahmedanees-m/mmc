"""Positive control for the calibration claim: does the edge-ablation gate discriminate?

The engineer-behavior result shows the loop's novel hypotheses conferred no held-out
predictive power. That alone does not show the gate is a discriminator rather than a blanket
rejecter: a gate that says no to everything is not demonstrated to say yes to anything true.
This is the positive control.

One strong-signal module holds the textbook Th2 circuit (GATA3 -> IL5 is the single strongest
knockdown effect in the atlas) alongside the demo's plausible-but-novel STK11 -> chemokine
hypothesis. The same edge-ablation gate (an edge is REQUIRED only if removing it lowers the
model's own held-out ACC_DEG, seed-stable) is run on every edge. Two clean outcomes:

  - Textbook edges are flagged required while the novel edges are not -> the gate
    discriminates true regulation from plausible-but-wrong, and the calibration claim is
    airtight.
  - Nothing is flagged required, even the textbook edges -> the gate is not shown to
    discriminate on this data; the claim is scoped to "the AI's plausible novel hypotheses
    conferred zero held-out predictive power" and drops "the gate supplies calibration."

Deterministic: the structure is hand-built (no model calls), only the store is queried.
Writes /app/gate_discrimination.json.
"""
from __future__ import annotations

import json

import numpy as np

from mmc.compile import structural
from mmc.eval.holdout import ModuleData, edge_ablation_holdout, loo_evaluate
from mmc.fit import fit_structural
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term
from mmc.shared import store

FDR = 0.10
COND = "Stim8hr"
SEEDS = (0, 1, 2)         # base held-out ACC_DEG is checked across these (fit is deterministic)
REQUIRED_DROP = 0.02      # a required edge lowers held-out ACC_DEG by at least this

GENES = ["STAT6", "GATA3", "TBX21", "STAT4", "IL4", "IL5", "IL13",
         "STK11", "CCL3", "CCL4", "CXCL8"]
TEXTBOOK = [("STAT6", "GATA3", 1), ("TBX21", "GATA3", -1), ("STAT4", "TBX21", 1),
            ("GATA3", "IL4", 1), ("GATA3", "IL5", 1), ("GATA3", "IL13", 1)]
NOVEL = [("STK11", "CCL3", -1), ("STK11", "CCL4", -1), ("STK11", "CXCL8", -1)]


def build_spec(genes, edges):
    E = [Edge(regulator=r, target=t, sign=s) for r, t, s in edges]
    by_target: dict[str, list] = {}
    for r, t, s in edges:
        by_target.setdefault(t, []).append((r, s))
    rules = {}
    for t, rs in by_target.items():
        names = [r for r, _ in rs][:3]                 # additive term, bounded to 3 regulators
        signs = {r: s for r, s in rs if r in names}
        rules[t] = Rule(terms=[Term(regulators=names, signs=signs)])
    return ModelSpec(genes=genes, edges=E, rules=rules)


def module_data(spec, condition):
    genes = list(spec.genes)
    gi = {g: i for i, g in enumerate(genes)}
    df = store.module_effects(genes, genes, condition)
    perts = sorted({p for p in df["perturbation"].tolist() if p in gi})
    pi = {p: i for i, p in enumerate(perts)}
    obs = np.zeros((len(perts), len(genes)))
    fdr = np.ones((len(perts), len(genes)))
    for _, r in df.iterrows():
        p, g = r["perturbation"], r["target_gene"]
        if p in pi and g in gi and p != g:
            obs[pi[p], gi[g]] = float(r["effect_size"])
            f = r["fdr"]
            fdr[pi[p], gi[g]] = float(f) if f == f else 1.0
    return ModuleData(genes, perts, obs, fdr < FDR, spec)


def backend():
    def fit_fn(spec, train_perts, train_obs, seed=0):
        observed = {p: {spec.genes[j]: float(train_obs[i, j])
                        for j in range(len(spec.genes)) if spec.genes[j] != p}
                    for i, p in enumerate(train_perts)}
        fits = fit_structural.multi_fit(spec, observed, n_starts=1, max_iter=150)
        return {"spec": spec, "params": fits[0]["params"]}

    def predict_fn(model, pert):
        return np.asarray(structural.knockdown(model["spec"], model["params"], pert))
    return fit_fn, predict_fn


def main():
    spec = build_spec(GENES, TEXTBOOK + NOVEL)
    mod = module_data(spec, COND)
    print(f"perturbations ({len(mod.perts)}): {mod.perts}", flush=True)
    fit_fn, predict_fn = backend()

    base = [loo_evaluate(mod, fit_fn, predict_fn, seed=s)["model"]["acc_deg"][0] for s in SEEDS]
    print(f"held-out ACC_DEG (model) by seed: {[round(b, 3) for b in base]}", flush=True)

    # edge-ablation drop per edge (seed 0; the fit is deterministic, base is seed-invariant)
    abl = {r["edge"]: r["drop"] for r in edge_ablation_holdout(mod, fit_fn, predict_fn, seed=0)}
    textbook_set = {f"{r}->{t}" for r, t, _ in TEXTBOOK}
    novel_set = {f"{r}->{t}" for r, t, _ in NOVEL}

    out = []
    for e in sorted(abl.keys()):
        drop = abl[e]
        cls = "textbook" if e in textbook_set else ("novel" if e in novel_set else "other")
        out.append({"edge": e, "class": cls, "drop": round(drop, 3),
                    "required": bool(drop >= REQUIRED_DROP)})
    out.sort(key=lambda r: -r["drop"])

    tb_req = [r["edge"] for r in out if r["class"] == "textbook" and r["required"]]
    nv_req = [r["edge"] for r in out if r["class"] == "novel" and r["required"]]
    result = {
        "perts": mod.perts,
        "held_out_acc_by_seed": [round(b, 3) for b in base],
        "required_drop_threshold": REQUIRED_DROP,
        "edges": out,
        "textbook_required": tb_req,
        "novel_required": nv_req,
        "gate_says_yes_to_a_true_edge": bool(tb_req),
        "gate_discriminates": bool(tb_req and not nv_req),
    }
    with open("/app/gate_discrimination.json", "w") as f:
        json.dump(result, f, indent=2)

    print("\n=== edge-ablation gate: held-out ACC_DEG drop when the edge is removed ===",
          flush=True)
    print(f"{'edge':<16}{'class':<10}{'held-out drop':<16}{'required'}")
    for r in out:
        print(f"{r['edge']:<16}{r['class']:<10}{r['drop']:<16}{r['required']}")
    print(f"\ntextbook edges flagged required: {tb_req}")
    print(f"novel edges flagged required:    {nv_req}")
    print(f"gate says yes to a true edge (positive control): {result['gate_says_yes_to_a_true_edge']}")
    print(f"gate discriminates true from plausible-but-novel: {result['gate_discriminates']}")
    print("wrote /app/gate_discrimination.json", flush=True)


if __name__ == "__main__":
    main()
