"""Scaled engineer-behavior characterization with confidence intervals.

Accumulates the loop's proposals across modules, conditions, and repeated runs, then measures:

  - validated rate: the fraction of novel edges that sit in a module whose model beats a
    linear baseline held-out;
  - modules-beating-linear rate: the fraction of module-condition models that beat linear;
  - reasoning versus search: the proposed structure against a random structure of equal edge
    count on in-sample fit (Pearson).

Rates are reported with Wilson score confidence intervals. Edge-level grounding is measured
separately in scripts/gate_discrimination.py. Results are written incrementally to
/app/engineer_behavior_scaled.json. Requires the model client (API key in the environment) and
the store.
"""
from __future__ import annotations

import json
import math

import numpy as np

from mmc.baselines import linear as linear_bl
from mmc.compile import structural
from mmc.data import module_extract
from mmc.eval.holdout import ModuleData, loo_evaluate
from mmc.fit import fit_structural
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term
from mmc.loop import run as runmod
from mmc.loop.run import discover
from mmc.shared import store

FDR = 0.10
OUT = "/app/engineer_behavior_scaled.json"
CYTO_CONTEXT = (
    "These genes are candidate regulators of cytokine production in stimulated CD4+ T cells "
    "together with the cytokine and chemokine outputs. Propose signed regulator-to-cytokine "
    "edges and a few regulator-to-regulator scaffold edges. Knockdown effects are log2 fold "
    "changes at 8 hours of stimulation."
)
GRID = [
    ("Cytokine_production", "Stim8hr", 3),
    ("Cytokine_production", "Stim48hr", 2),
    ("Th2_GATA3", "Stim8hr", 2),
    ("CD4_TF_network", "Stim8hr", 1),
    ("TCR_core", "Stim8hr", 1),
]


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


def backend(n_starts=1):
    def fit_fn(spec, train_perts, train_obs, seed=0):
        observed = {p: {spec.genes[j]: float(train_obs[i, j])
                        for j in range(len(spec.genes)) if spec.genes[j] != p}
                    for i, p in enumerate(train_perts)}
        fits = fit_structural.multi_fit(spec, observed, n_starts=n_starts, max_iter=120)
        return {"spec": spec, "params": fits[0]["params"]}

    def predict_fn(model, pert):
        return np.asarray(structural.knockdown(model["spec"], model["params"], pert))
    return fit_fn, predict_fn


def linear_fn_for(genes):
    def linear_fn(train_perts, train_obs, held):
        td = {p: {genes[j]: float(train_obs[i, j]) for j in range(len(genes))}
              for i, p in enumerate(train_perts)}
        pred = linear_bl.reconstruct(td, genes, list(train_perts), [held])[held]
        return np.array([pred.get(g, 0.0) for g in genes])
    return linear_fn


def wilson(k, n):
    if n == 0:
        return [None, None, None]
    p = k / n
    z = 1.96
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(p, 3), round(max(0.0, centre - half), 3), round(min(1.0, centre + half), 3)]


def random_structure(genes, n_edges, rng):
    edges, by_t, seen = [], {}, set()
    tries = 0
    while len(edges) < n_edges and tries < n_edges * 30 + 50:
        tries += 1
        r, t = str(rng.choice(genes)), str(rng.choice(genes))
        if r == t or (r, t) in seen or len(by_t.get(t, [])) >= 9:
            continue
        s = int(rng.choice([1, -1]))
        edges.append(Edge(regulator=r, target=t, sign=s))
        seen.add((r, t))
        by_t.setdefault(t, []).append((r, s))
    rules = {}
    for t, rs in by_t.items():
        terms = []
        for k in range(0, len(rs), 3):
            ch = rs[k:k + 3]
            terms.append(Term(regulators=[r for r, _ in ch], signs={r: sg for r, sg in ch}))
            if len(terms) == 3:
                break
        rules[t] = Rule(terms=terms)
    return ModelSpec(genes=list(genes), edges=edges, rules=rules)


def insample_pearson(spec, condition):
    mod = module_data(spec, condition)
    fit_fn, predict_fn = backend(n_starts=2)
    model = fit_fn(spec, mod.perts, mod.observed)
    pred = np.array([predict_fn(model, p) for p in mod.perts])
    a, b = pred.ravel(), mod.observed.ravel()
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def main():
    with open("/app/cytokine_module_def.json") as f:
        cyto = json.load(f)
    module_extract.register_module("Cytokine_production", cyto["regulators"], cyto["targets"])
    runmod.DEFAULT_CONTEXT["Cytokine_production"] = CYTO_CONTEXT
    dark_for = {"Cytokine_production": set(cyto["dark"])}

    result = {"grid": [list(g) for g in GRID], "n_proposals": 0,
              "novel_edges": [], "modules": []}
    rng = np.random.default_rng(0)

    for name, condition, n_runs in GRID:
        dark = dark_for.get(name, set())
        for run_i in range(n_runs):
            try:
                dis = discover(name, condition, backend="structural", max_iters=3,
                               n_starts=2, max_iter=120)
            except Exception as e:  # noqa: BLE001
                result["modules"].append({"module": name, "condition": condition, "run": run_i,
                                          "error": str(e)[:200]})
                continue
            props = [s.rationale for s in dis.trace.steps
                     if s.kind in ("propose", "repair") and s.rationale]
            result["n_proposals"] += len(props)
            frozen = dis.ensemble.best().spec
            novel = [(e.regulator, e.target, e.sign) for e in frozen.edges
                     if e.regulator in dark]

            mod = module_data(frozen, condition)
            fit_fn, predict_fn = backend(n_starts=1)
            rep = loo_evaluate(mod, fit_fn, predict_fn, linear_fn=linear_fn_for(mod.genes),
                               seed=0)
            mde = rep.get("model", {}).get("de_overlap")
            lde = rep.get("linear", {}).get("de_overlap")
            beats = bool(mde and lde and mde[0] is not None and lde[0] is not None
                         and mde[0] > lde[0])
            for r, t, _s in novel:
                result["novel_edges"].append(
                    {"module": name, "condition": condition, "run": run_i, "edge": f"{r}->{t}",
                     "in_beats_linear_module": beats})

            try:
                rand = random_structure(list(frozen.genes), len(frozen.edges), rng)
                claude_p = insample_pearson(frozen, condition)
                rand_p = insample_pearson(rand, condition)
            except Exception:  # noqa: BLE001
                claude_p = rand_p = None

            result["modules"].append(
                {"module": name, "condition": condition, "run": run_i,
                 "n_proposals": len(props), "n_novel": len(novel), "beats_linear": beats,
                 "model_de": mde, "linear_de": lde,
                 "claude_insample_pearson": claude_p, "random_insample_pearson": rand_p})
            with open(OUT, "w") as f:
                json.dump(result, f, indent=2)
            print(f"{name}/{condition}/run{run_i}: props={len(props)} novel={len(novel)} "
                  f"beats_linear={beats} claudeP={claude_p} randP={rand_p}", flush=True)

    ne = result["novel_edges"]
    validated = sum(1 for e in ne if e["in_beats_linear_module"])
    mods = [m for m in result["modules"] if "error" not in m]
    beat = sum(1 for m in mods if m["beats_linear"])
    cp = [m["claude_insample_pearson"] for m in mods
          if m.get("claude_insample_pearson") is not None]
    rp = [m["random_insample_pearson"] for m in mods
          if m.get("random_insample_pearson") is not None]
    result["summary"] = {
        "n_runs": len(mods),
        "n_proposals": result["n_proposals"],
        "n_novel_edges": len(ne),
        "validated_rate": wilson(validated, len(ne)),
        "modules_beating_linear": [beat, len(mods)],
        "beats_linear_rate": wilson(beat, len(mods)),
        "claude_insample_pearson_mean": round(float(np.mean(cp)), 3) if cp else None,
        "random_insample_pearson_mean": round(float(np.mean(rp)), 3) if rp else None,
    }
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)
    print("\n=== SUMMARY ===")
    print(json.dumps(result["summary"], indent=2), flush=True)


if __name__ == "__main__":
    main()
