"""Pre-compute the interactive-circuit responses for the demo.

Builds the ablation-validated Th2 / GATA3 sub-circuit (the edges the gate-discrimination
control flagged as held-out-necessary), fits it to the store, and simulates the structural
model's response to a small set of interventions restricted to validated regulators, including
a double knockdown. For each intervention it records the per-gene response, the causal path
(the reachable validated edges), and the additive baseline (the observed single knockdowns, and
their sum for the double). The demo reads the resulting JSON and needs no store or backend.

Writes /app/demo_interventions.json. Requires the store; deterministic, no model calls.
"""
from __future__ import annotations

import json

import numpy as np

from mmc.compile import structural
from mmc.fit import fit_structural
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term
from mmc.shared import store

COND = "Stim8hr"
# validated Th2 sub-circuit (gate_discrimination held-out-necessary edges)
GENES = ["STAT4", "TBX21", "STAT6", "GATA3", "IL4", "IL5", "IL13"]
EDGES = [("STAT6", "GATA3", 1), ("TBX21", "GATA3", -1), ("STAT4", "TBX21", 1),
         ("GATA3", "IL4", 1), ("GATA3", "IL5", 1), ("GATA3", "IL13", 1)]
READOUTS = ["IL4", "IL5", "IL13"]
INTERVENTIONS = [
    ("GATA3 knockdown", "knockdown", ["GATA3"]),
    ("GATA3 activation", "activation", ["GATA3"]),
    ("STAT6 knockdown", "knockdown", ["STAT6"]),
    ("GATA3 + STAT6 knockdown", "knockdown", ["GATA3", "STAT6"]),
]


def build_spec():
    edges = [Edge(regulator=r, target=t, sign=s) for r, t, s in EDGES]
    by_t = {}
    for r, t, s in EDGES:
        by_t.setdefault(t, []).append((r, s))
    rules = {t: Rule(terms=[Term(regulators=[r for r, _ in rs],
                                 signs={r: s for r, s in rs})])
             for t, rs in by_t.items()}
    return ModelSpec(genes=list(GENES), edges=edges, rules=rules)


def observed_single(gene):
    """Observed knockdown delta of `gene` over the circuit genes, from the store."""
    df = store.module_effects([gene], GENES, COND)
    out = {g: 0.0 for g in GENES}
    for _, r in df.iterrows():
        if r["perturbation"] == gene and r["target_gene"] in out:
            out[r["target_gene"]] = float(r["effect_size"])
    return out


def causal_path(sources):
    """Validated edges reachable from the intervened genes (BFS over the DAG)."""
    reach, frontier, used = set(sources), list(sources), []
    while frontier:
        g = frontier.pop()
        for r, t, _s in EDGES:
            if r == g and (r, t) not in used:
                used.append((r, t))
                if t not in reach:
                    reach.add(t)
                    frontier.append(t)
    return [[r, t] for r, t in used]


def main():
    spec = build_spec()
    genes = list(spec.genes)
    gi = {g: i for i, g in enumerate(genes)}
    df = store.module_effects(genes, genes, COND)
    perts = sorted({p for p in df["perturbation"].tolist() if p in gi})
    observed = {}
    for _, r in df.iterrows():
        p, g = r["perturbation"], r["target_gene"]
        if p in perts and g in gi and p != g:
            observed.setdefault(p, {})[g] = float(r["effect_size"])
    fits = fit_structural.multi_fit(spec, observed, n_starts=3, max_iter=200)
    params = fits[0]["params"]

    obs_singles = {g: observed_single(g) for g in ("GATA3", "STAT6")}

    out = {"genes": genes, "edges": [[r, t, s] for r, t, s in EDGES],
           "readouts": READOUTS, "interventions": {}}
    for name, kind, targets in INTERVENTIONS:
        if len(targets) == 1 and kind == "knockdown":
            delta = np.asarray(structural.knockdown(spec, params, targets[0]))
        elif len(targets) == 1 and kind == "activation":
            delta = np.asarray(structural.activation(spec, params, targets[0]))
        else:
            level = 0.0 if kind == "knockdown" else structural.ACTIVATION_LEVEL
            delta = np.asarray(structural.perturb_set(spec, params, targets, level))
        response = {g: round(float(delta[gi[g]]), 3) for g in genes}
        if kind == "knockdown":
            add = {g: round(sum(obs_singles[t][g] for t in targets), 3) for g in genes}
        else:
            add = None
        out["interventions"][name] = {
            "kind": kind, "targets": targets,
            "response": response, "path": causal_path(targets),
            "additive_baseline": add,
            "composed": len(targets) > 1}

    with open("/app/demo_interventions.json", "w") as f:
        json.dump(out, f, indent=2)
    print("wrote /app/demo_interventions.json")
    for name, d in out["interventions"].items():
        ro = {g: d["response"][g] for g in READOUTS}
        print(f"  {name}: readouts {ro}  path {len(d['path'])} edges")


if __name__ == "__main__":
    main()
