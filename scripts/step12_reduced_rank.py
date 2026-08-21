"""Amendment A18: does a reduced-rank map recover the full ridge map's held-out score?

Claim 2 rests on a property of the response matrix measured by singular value
decomposition: an effective rank of 3.64 out of 28 perturbations on the primary module,
which no sparse causal generator reproduces. That is a statement about the data, and a
diagnostic is easier to argue with than a prediction.

This script converts it into one. If the response matrix is effectively low rank, then
truncating the ridge map near that rank should cost almost nothing on held-out folds.
The prediction was fixed in the pre-registration before this ran: across the 13 modules
where the linear arm beats its own random null, the full ridge clears the reduced-rank
map on at most 4 modules by the section 1.2 statistic, and the median ratio of
reduced-rank to full-rank held-out DE-overlap is at least 0.90. Refuted if the ridge
clears on 5 or more, or if the median ratio falls below 0.90.

The rank is chosen per fold from the training rows alone, so no fold's own response
informs the rank used to predict it. Fixed ranks of 3 and 4 are reported alongside as
pre-specified alternatives, and the sweep from 1 to 8 is descriptive only.

Per-fold scores are written out for every arm, which the earlier comparator runs did not
do (section 2.31).

    python scripts/step12_reduced_rank.py --module coresponse_PIM1 \\
        --module-def /work/step7_defs/coresponse_PIM1.json \\
        --out /work/results/step12/coresponse_PIM1.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mmc.data import module_data, module_extract
from mmc.eval import compare
from mmc.eval import identifiability as ident

# Descriptive only, and explicitly outside the A18 test. Extended beyond the 1 to 8 named
# in the amendment on 2026-08-22, after the first module returned an effective rank of 8.1
# out of 40 perturbations: a curve that stops at the adaptive rank cannot show where it
# converges. The pre-registered rule is untouched, since it is stated on the adaptive arm
# and on fixed ranks 3 and 4.
SWEEP_BASE = (1, 2, 3, 4, 5, 6, 7, 8, 12, 16, 24, 32)


def arm(mod, name: str, rank) -> dict:
    """Score one linear arm over all leave-one-out folds."""
    fn = (compare.bind_linear(mod) if rank == "full"
          else compare.bind_reduced_rank(mod, rank=rank))
    preds = compare.fold_predictions(mod, fn)
    scores = compare.score_predictions(mod, preds)
    return {"arm": name, "rank": rank,
            "de_overlap": compare.bootstrap_mean(scores["de_overlap"]),
            "acc_deg": compare.bootstrap_mean(scores["acc_deg"]),
            "per_fold": {"de_overlap": [None if np.isnan(v) else float(v)
                                        for v in scores["de_overlap"]],
                         "acc_deg": [None if np.isnan(v) else float(v)
                                     for v in scores["acc_deg"]]},
            "_folds": scores["de_overlap"]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--module-def", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(d["module"], d["regulators"], d["targets"])

    mod = module_data.build_module_data(args.module, args.condition)
    n = len(mod.perts)

    # the per-fold ranks the adaptive arm will use, recorded so the choice is inspectable
    ranks = [int(np.ceil(ident.effective_rank(mod.observed[[j for j in range(n) if j != i]])))
             for i in range(n)]
    full_rank = float(ident.effective_rank(mod.observed))
    print(f"=== {args.module}: {len(mod.genes)} genes, {n} perturbations ===")
    print(f"effective rank on the full matrix {full_rank:.3f}; "
          f"per-fold adaptive rank {min(ranks)} to {max(ranks)}", flush=True)

    arms = {"full": arm(mod, "ridge, full rank", "full"),
            "adaptive": arm(mod, "reduced rank, per-fold", None)}
    for k in [k for k in SWEEP_BASE if k < n]:
        arms[f"rank_{k}"] = arm(mod, f"reduced rank {k}", k)

    ridge = arms["full"].pop("_folds")
    out = {"module": args.module, "condition": args.condition,
           "n_perturbations": n, "n_genes": len(mod.genes),
           "effective_rank_full": full_rank,
           "adaptive_ranks": ranks,
           "arms": {}, "tests": {}}

    for key, a in arms.items():
        folds = a.pop("_folds", None)
        out["arms"][key] = a
        if folds is not None:
            # does the full ridge clear this reduced-rank arm? advantage True means it does,
            # which counts against the A18 prediction
            d = compare.paired_delta(ridge, folds)
            ratio = (a["de_overlap"]["mean"] / arms["full"]["de_overlap"]["mean"]
                     if arms["full"]["de_overlap"]["mean"] else float("nan"))
            out["tests"][key] = {"ridge_minus_reduced": d,
                                 "ridge_clears": bool(d["advantage"]),
                                 "ratio_reduced_over_full": float(ratio)}

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)

    r = out["tests"]["adaptive"]
    print(f"ridge {arms['full']['de_overlap']['mean']:.4f}  "
          f"reduced (adaptive) {arms['adaptive']['de_overlap']['mean']:.4f}  "
          f"ratio {r['ratio_reduced_over_full']:.3f}")
    print(f"ridge minus reduced {r['ridge_minus_reduced']['delta']:+.4f} "
          f"[{r['ridge_minus_reduced']['lo']:+.4f}, {r['ridge_minus_reduced']['hi']:+.4f}]  "
          f"ridge clears: {r['ridge_clears']}")
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
