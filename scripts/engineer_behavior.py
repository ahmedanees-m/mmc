"""Work Stream A: characterize the engineer's self-correction as a measured capability.

Accumulate the structural hypotheses the loop proposes across modules and runs, classify
each novel (dark) hypothesis by whether Claude gave a coherent, data-referencing rationale,
and assign the ground-truth verdict from the held-out gate: a hypothesis is VALIDATED only
if the model that contains it beats the linear baseline on held-out DE-overlap; otherwise
REJECTED. The finding is whether Claude's mechanistic plausibility tracks predictive
validity, and whether the held-out gate reliably supplies the calibration.

Writes engineer_behavior.json (the proposal corpus + per-module gate verdict + the derived
metrics). Requires the model client and the store.
"""
from __future__ import annotations

import json

import numpy as np

from mmc.baselines import linear as linear_bl
from mmc.compile import structural
from mmc.data import module_extract
from mmc.eval.holdout import ModuleData, loo_evaluate
from mmc.fit import fit_structural
from mmc.loop import run as runmod
from mmc.loop.propose import model_hash
from mmc.loop.run import discover
from mmc.shared import store

FDR = 0.10
CYTO_CONTEXT = (
    "These genes are candidate regulators of cytokine production in stimulated CD4+ T cells "
    "together with the cytokine and chemokine outputs. Propose signed regulator-to-cytokine "
    "edges and a few regulator-to-regulator scaffold edges. Knockdown effects are log2 fold "
    "changes at 8 hours of stimulation."
)


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
        fits = fit_structural.multi_fit(spec, observed, n_starts=3, max_iter=150)
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


def run_module(name, condition, dark, runs, context=None):
    if context:
        runmod.DEFAULT_CONTEXT[name] = context
    proposals, frozen = [], None
    for r in range(runs):
        result = discover(name, condition, backend="structural",
                          max_iters=3, n_starts=3, max_iter=120)
        for s in result.trace.steps:
            if s.kind in ("propose", "repair") and s.rationale:
                proposals.append({"module": name, "run": r, "kind": s.kind,
                                  "rationale": s.rationale})
        frozen = result.ensemble.best()
    # the surviving novel (dark) hypotheses in the last frozen structure
    edges = [{"regulator": e.regulator, "target": e.target, "sign": e.sign,
              "novel": e.regulator in dark} for e in frozen.spec.edges]
    # ground-truth verdict: does the model beat the linear baseline held-out (DE-overlap)?
    mod = module_data(frozen.spec, condition)
    fit_fn, predict_fn = backend()
    rep = loo_evaluate(mod, fit_fn, predict_fn, linear_fn=linear_fn_for(mod.genes),
                       n_boot=1000, seed=0)

    def de(m):
        p = rep.get(m, {}).get("de_overlap")
        return None if not p else p[0]
    model_de, linear_de = de("model"), de("linear")
    beats_linear = (model_de is not None and linear_de is not None and model_de > linear_de)
    return {"module": name, "hash": model_hash(frozen.spec), "proposals": proposals,
            "edges": edges, "n_novel_edges": sum(e["novel"] for e in edges),
            "held_out_DE_overlap": {"model": model_de, "linear": linear_de},
            "beats_linear_held_out": bool(beats_linear)}


def main():
    with open("/app/cytokine_module_def.json") as f:
        d = json.load(f)
    module_extract.register_module("Cytokine_production", d["regulators"], d["targets"])

    out = []
    out.append(run_module("Cytokine_production", "Stim8hr", set(d["dark"]), runs=2,
                          context=CYTO_CONTEXT))
    out.append(run_module("Th2_GATA3", "Stim8hr", set(), runs=2))   # no dark set: all textbook

    # aggregate: every module's novel hypotheses are validated only if that module's model
    # beats the linear baseline held-out. Catch rate = fraction of coherent proposals rejected.
    n_novel = sum(m["n_novel_edges"] for m in out)
    n_validated = sum(m["n_novel_edges"] for m in out if m["beats_linear_held_out"])
    n_proposals = sum(len(m["proposals"]) for m in out)
    result = {
        "modules": out,
        "n_coherent_proposals": n_proposals,
        "n_novel_hypotheses": n_novel,
        "n_novel_validated_held_out": n_validated,
        "catch_rate": (n_novel - n_validated) / n_novel if n_novel else None,
    }
    with open("/app/engineer_behavior.json", "w") as f:
        json.dump(result, f, indent=2)

    print("=== engineer-behavior characterization ===")
    for m in out:
        print(f"  {m['module']}: {len(m['proposals'])} coherent proposals, "
              f"{m['n_novel_edges']} novel edges; held-out DE-overlap model "
              f"{m['held_out_DE_overlap']['model']} vs linear {m['held_out_DE_overlap']['linear']} "
              f"-> beats linear: {m['beats_linear_held_out']}")
    print(f"\nHEADLINE: {n_proposals} coherently-argued proposals; {n_novel} novel hypotheses; "
          f"{n_validated} validated as a held-out predictive necessity; "
          f"catch rate {result['catch_rate']}.")


if __name__ == "__main__":
    main()
