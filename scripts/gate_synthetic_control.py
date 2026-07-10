"""Synthetic positive control: the module-level held-out gate is a verified discriminator.

Data is generated from a known sparse true structure with additive noise. Leave-one-
perturbation-out, we fit the TRUE structure and a fully-connected (over-connected) structure
and score each against the mean baseline on held-out DE-overlap. The claim to test: the true
structure clears the held-out-advantage bar and the over-connected one does not, so the
module-level gate can fire positive and rejects over-fit structure, meaning its rejection of
the real modules is discrimination, not a uniform no.

Deterministic (seeded noise); no model calls. Writes /app/gate_synthetic_control.json.
"""
from __future__ import annotations

import json

import numpy as np

from mmc.compile.structural import knockdown, steady_state
from mmc.eval.holdout import ModuleData, loo_evaluate
from mmc.fit import fit_structural
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term

GENES = ["A", "B", "C", "D", "E"]
TRUE_EDGES = [("A", "B", 1), ("B", "C", -1), ("A", "C", 1), ("C", "D", 1), ("D", "E", -1)]
# fully-connected forward DAG in the A<B<C<D<E order: maximally over-connected, no loops
FULL_EDGES = [(GENES[i], GENES[j], 1) for j in range(len(GENES)) for i in range(j)]
NOISE = 0.2
DE_T = 0.15

TRUE_PARAMS = {
    "basal": np.array([1.0, 0.2, 0.2, 0.2, 0.2]),
    "terms": {"B": [{"prod": 2.0, "w": {"A": 4.0}, "theta": {"A": 1.0}}],
              "C": [{"prod": 2.0, "w": {"B": 4.0}, "theta": {"B": 1.0}},
                    {"prod": 1.5, "w": {"A": 3.0}, "theta": {"A": 1.0}}],
              "D": [{"prod": 2.0, "w": {"C": 4.0}, "theta": {"C": 1.0}}],
              "E": [{"prod": 2.0, "w": {"D": 4.0}, "theta": {"D": 1.0}}]},
}


def build_spec(genes, edges):
    E = [Edge(regulator=r, target=t, sign=s) for r, t, s in edges]
    by_target: dict[str, list] = {}
    for r, t, s in edges:
        by_target.setdefault(t, []).append((r, s))
    rules = {}
    for t, rs in by_target.items():
        rs = rs[:9]
        terms = []
        for k in range(0, len(rs), 3):
            chunk = rs[k:k + 3]
            terms.append(Term(regulators=[r for r, _ in chunk],
                              signs={r: s for r, s in chunk}))
            if len(terms) == 3:
                break
        rules[t] = Rule(terms=terms)
    return ModelSpec(genes=genes, edges=E, rules=rules)


def backend():
    def fit_fn(spec, tp, tobs, seed=0):
        obs = {p: {spec.genes[j]: float(tobs[i, j])
                   for j in range(len(spec.genes)) if spec.genes[j] != p}
               for i, p in enumerate(tp)}
        fits = fit_structural.multi_fit(spec, obs, n_starts=4, max_iter=300)
        return {"spec": spec, "params": fits[0]["params"]}

    def predict_fn(m, p):
        return np.asarray(knockdown(m["spec"], m["params"], p))
    return fit_fn, predict_fn


def main():
    true_spec = build_spec(GENES, TRUE_EDGES)
    wt = steady_state(true_spec, TRUE_PARAMS)
    clean = np.zeros((len(GENES), len(GENES)))
    for i, reg in enumerate(GENES):
        clean[i] = knockdown(true_spec, TRUE_PARAMS, reg, wt=wt)
    rng = np.random.default_rng(0)
    obs = clean + rng.normal(0, NOISE, clean.shape)
    de = np.abs(clean) > DE_T                       # DE mask from the clean ground truth
    print(f"DE genes per perturbation: {de.sum(axis=1).tolist()}", flush=True)
    fit_fn, predict_fn = backend()

    out = {}
    for name, edges in [("true_sparse", TRUE_EDGES), ("over_connected_full", FULL_EDGES)]:
        spec = build_spec(GENES, edges)
        mod = ModuleData(GENES, GENES, obs, de, spec)
        rep = loo_evaluate(mod, fit_fn, predict_fn, seed=0)
        md, mn = rep["model"]["de_overlap"], rep["mean"]["de_overlap"]
        ma, mna = rep["model"]["acc_deg"], rep["mean"]["acc_deg"]
        beats = md[0] > mn[0]
        out[name] = {"n_edges": len(edges),
                     "model_de_overlap": [round(x, 3) for x in md],
                     "mean_de_overlap": [round(x, 3) for x in mn],
                     "model_acc_deg": [round(x, 3) for x in ma],
                     "mean_acc_deg": [round(x, 3) for x in mna],
                     "beats_mean_de_overlap": bool(beats)}
        print(f"{name} ({len(edges)} edges): model DE-overlap {md[0]:.3f} "
              f"[{md[1]:.3f},{md[2]:.3f}]  mean {mn[0]:.3f}  beats_mean {beats}", flush=True)

    t, o = out["true_sparse"], out["over_connected_full"]
    true_certifies = t["model_de_overlap"][1] > t["mean_de_overlap"][0]     # true CI-low > mean point
    over_below_true = o["model_de_overlap"][0] < t["model_de_overlap"][0] - 0.15
    out["true_certifies_separated"] = bool(true_certifies)
    out["over_connected_ranked_below_true"] = bool(over_below_true)
    out["gate_discriminates"] = bool(true_certifies and over_below_true)
    with open("/app/gate_synthetic_control.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\ntrue certifies (CI-low above mean): {true_certifies}; over-connected ranked below "
          f"true: {over_below_true}; gate discriminates: {out['gate_discriminates']}", flush=True)


if __name__ == "__main__":
    main()
