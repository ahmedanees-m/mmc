"""Step 3, calibration gate: can a simulator reproduce the real module's identifiability?

PREREG_v4 section 4.1 makes the power curve conditional on a stated acceptance
criterion: the simulated module must match the real one on effective rank within 15
percent relative and on leading-PC variance fraction within 0.05 absolute. Both, not
either. If that fails, section 4.1 says the mismatch is itself reported and the phase
diagram is scoped down to a qualitative boundary.

The Step 1 decomposition predicts this will be tight. Seventy-eight percent of what the
linear baseline achieves on the cytokine module comes from the mean training response,
so the real response matrix is dominated by a shared component, and its effective rank
is 3.64 out of 28 perturbations. A generator that simulates each knockdown through a
sparse causal graph produces responses supported only on each perturbed gene's
descendants, which are close to independent across perturbations and therefore
higher-rank. The question this script answers is how much of a shared component has to
be added before the simulator matches, and whether the structural part survives it.

Three generators are compared against the real module:

  structural        knockdowns simulated through a fitted structural model, plus
                    bootstrapped residual noise
  structural+shared the same, plus a shared response direction at a swept weight
  lowrank           a pure low-rank plus noise generator with no causal structure,
                    as the null model for the diagnostics themselves

    python scripts/step3_calibrate.py --module Cytokine_production \\
        --module-def /work/cytokine_module_def.json --out /work/results/step3_calib.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mmc.compile import structural
from mmc.data import module_data, module_extract
from mmc.eval import identifiability as ident
from mmc.eval import random_null
from mmc.fit import fit_structural

RANK_TOL_REL = 0.15        # section 4.1
PC1_TOL_ABS = 0.05


def diagnostics_of(matrix: np.ndarray) -> dict:
    return {
        "effective_rank": ident.effective_rank(matrix),
        "leading_pc_fraction": ident.mean_leading_pc_fraction(matrix),
        "perturbation_specific_ratio": ident.perturbation_specific_ratio(matrix),
    }


def accepts(sim: dict, real: dict) -> dict:
    rank_rel = abs(sim["effective_rank"] - real["effective_rank"]) / real["effective_rank"]
    pc1_abs = abs(sim["leading_pc_fraction"] - real["leading_pc_fraction"])
    return {
        "effective_rank_relative_error": float(rank_rel),
        "leading_pc_absolute_error": float(pc1_abs),
        "rank_ok": bool(rank_rel <= RANK_TOL_REL),
        "pc1_ok": bool(pc1_abs <= PC1_TOL_ABS),
        "accepted": bool(rank_rel <= RANK_TOL_REL and pc1_abs <= PC1_TOL_ABS),
    }


def fit_generator(mod, n_edges: int, seed: int):
    """A structural model fitted to the real module, used as the data generator."""
    genes = list(mod.genes)
    rng = np.random.default_rng(seed)
    spec = random_null.sample_spec(genes, mod.perts, n_edges, rng)
    observed = {
        mod.perts[i]: {genes[j]: float(mod.observed[i, j])
                       for j in range(len(genes)) if genes[j] != mod.perts[i]}
        for i in range(len(mod.perts))
    }
    fits = fit_structural.multi_fit(spec, observed, n_starts=4, max_iter=250)
    return spec, fits[0]["params"]


def simulate(mod, spec, params, *, noise_scale: float, shared_weight: float,
             seed: int, shared_rank: int = 2) -> np.ndarray:
    """Simulate the module's response matrix from the fitted structural model.

    `shared_weight` mixes in a low-rank co-response component taken from the real
    module. It has to be a *per-perturbation scaled* component, not a constant vector
    added to every row: both diagnostics are computed on the column-centred matrix, so
    a constant offset is removed by the centring and cannot change either of them. An
    earlier version of this script added a constant and produced identical diagnostics
    at every weight, which is the signature of exactly that mistake.
    """
    rng = np.random.default_rng(seed)
    clean = np.stack([np.asarray(structural.knockdown(spec, params, p))
                      for p in mod.perts])

    # residual noise bootstrapped from the real module rather than assumed Gaussian
    resid = np.asarray(mod.observed, float) - clean
    noise = rng.choice(resid.reshape(-1), size=clean.shape, replace=True) * noise_scale

    out = clean + noise
    if shared_weight > 0:
        real = np.asarray(mod.observed, float)
        centred = real - real.mean(axis=0, keepdims=True)
        u, s, vt = np.linalg.svd(centred, full_matrices=False)
        k = min(shared_rank, s.size)
        out = out + shared_weight * ((u[:, :k] * s[:k]) @ vt[:k])
    return out


def lowrank_surrogate(mod, rank: int, seed: int) -> np.ndarray:
    """Low rank plus noise, with no causal structure. The null model for the diagnostics."""
    rng = np.random.default_rng(seed)
    real = np.asarray(mod.observed, float)
    u, s, vt = np.linalg.svd(real - real.mean(axis=0, keepdims=True), full_matrices=False)
    approx = (u[:, :rank] * s[:rank]) @ vt[:rank]
    resid = real - real.mean(axis=0, keepdims=True) - approx
    noise = rng.choice(resid.reshape(-1), size=real.shape, replace=True)
    return real.mean(axis=0, keepdims=True) + approx + noise


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", default="Cytokine_production")
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--module-def", default="")
    ap.add_argument("--edges", type=int, default=24)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(d["module"], d["regulators"], d["targets"])

    mod = module_data.build_module_data(args.module, args.condition)
    real = diagnostics_of(np.asarray(mod.observed, float))
    print(f"=== {args.module}: {len(mod.genes)} genes, {len(mod.perts)} perturbations ===")
    print(f"real: effective rank {real['effective_rank']:.3f}, "
          f"leading-PC {real['leading_pc_fraction']:.4f}, "
          f"specific/shared {real['perturbation_specific_ratio']:.3f}")
    print(f"acceptance: rank within {RANK_TOL_REL:.0%} relative, "
          f"leading-PC within {PC1_TOL_ABS} absolute\n")

    out = {"module": args.module, "condition": args.condition, "real": real,
           "tolerances": {"rank_relative": RANK_TOL_REL, "pc1_absolute": PC1_TOL_ABS},
           "variants": []}

    print("fitting the structural generator to the real module")
    spec, params = fit_generator(mod, args.edges, seed=0)
    print(f"  generator: {len(spec.edges)} edges\n")

    grid = [("structural", 0.0), ("structural+lowrank", 0.5), ("structural+lowrank", 1.0),
            ("structural+lowrank", 2.0), ("structural+lowrank", 4.0),
            ("structural+lowrank", 8.0)]
    for name, w in grid:
        sims = [diagnostics_of(simulate(mod, spec, params, noise_scale=1.0,
                                        shared_weight=w, seed=s))
                for s in range(args.seeds)]
        avg = {k: float(np.mean([s[k] for s in sims])) for k in sims[0]}
        verdict = accepts(avg, real)
        out["variants"].append({"generator": name, "shared_weight": w,
                                "diagnostics": avg, **verdict})
        print(f"  {name:<20} shared={w:<5} rank {avg['effective_rank']:6.3f} "
              f"(err {verdict['effective_rank_relative_error']:5.1%})  "
              f"PC1 {avg['leading_pc_fraction']:.4f} "
              f"(err {verdict['leading_pc_absolute_error']:.4f})  "
              f"{'ACCEPT' if verdict['accepted'] else 'reject'}")

    for rank in (1, 2, 3, 4, 6):
        sims = [diagnostics_of(lowrank_surrogate(mod, rank, s)) for s in range(args.seeds)]
        avg = {k: float(np.mean([s[k] for s in sims])) for k in sims[0]}
        verdict = accepts(avg, real)
        out["variants"].append({"generator": "lowrank", "rank": rank,
                                "diagnostics": avg, **verdict})
        print(f"  {'lowrank':<20} rank={rank:<5} rank {avg['effective_rank']:6.3f} "
              f"(err {verdict['effective_rank_relative_error']:5.1%})  "
              f"PC1 {avg['leading_pc_fraction']:.4f} "
              f"(err {verdict['leading_pc_absolute_error']:.4f})  "
              f"{'ACCEPT' if verdict['accepted'] else 'reject'}")

    accepted = [v for v in out["variants"] if v["accepted"]]
    out["any_accepted"] = bool(accepted)
    out["accepted_variants"] = [{k: v[k] for k in ("generator", "shared_weight", "rank")
                                 if k in v} for v in accepted]
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2, default=float)

    print()
    if accepted:
        print(f"{len(accepted)} variant(s) meet the section 4.1 criterion; the power "
              f"curve proceeds quantitatively on those.")
    else:
        print("No variant meets the section 4.1 criterion. Per section 4.1 the power "
              "curve is scoped down to a qualitative boundary and the mismatch is "
              "reported as a finding about current simulators.")
    print(f"wrote {outpath}")


if __name__ == "__main__":
    main()
