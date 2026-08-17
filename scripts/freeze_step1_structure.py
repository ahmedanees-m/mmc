"""Freeze an S1 structure for the Step 1 comparator (PREREG_v4 section 2.2).

Runs the proposal loop on a module and writes the frozen structure in the shape
`step1_comparator.py --frozen` reads. Kept separate from `freeze_model.py`, which
freezes against the superseded two-tier transfer split; Step 1 evaluates
leave-one-perturbation-out over all of a module's perturbations, so the structure is
proposed from the module as a whole.

That the proposal sees the whole module matters for how the S1 row is read, and it is
recorded in the output rather than left implicit. The harness refits parameters inside
every fold, but the structure is fixed across folds for every source, so S1 has seen
the held-out perturbations at proposal time in the same sense that the mean-difference
and GRNBoost2 arms have. The textbook arm is the only one whose edges predate this
data entirely, and the nested oracle is the only arm with an outer split that
selection never touched.

    python scripts/freeze_step1_structure.py Cytokine_production Stim8hr \\
        --module-def /work/cytokine_module_def.json --out /work/frozen_cytokine.json

Requires MMC_ZHU_STORE and ANTHROPIC_API_KEY. Every call is added to the ledger in
paper/LEDGER_api_spend.md by hand after the run.
"""
from __future__ import annotations

import argparse
import json

from mmc.data import module_extract
from mmc.loop.run import discover


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("module")
    ap.add_argument("condition", nargs="?", default="Stim8hr")
    ap.add_argument("--module-def", default="")
    ap.add_argument("--max-iters", type=int, default=3)
    ap.add_argument("--n-starts", type=int, default=6)
    ap.add_argument("--max-iter", type=int, default=250)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.module_def:
        with open(args.module_def) as f:
            d = json.load(f)
        module_extract.register_module(args.module, d["regulators"], d["targets"])

    result = discover(args.module, args.condition, max_iters=args.max_iters,
                      n_starts=args.n_starts, max_iter=args.max_iter)
    best = result.ensemble.best()
    spec = best.spec

    payload = {
        "module": args.module,
        "state": args.condition,
        "genes": list(spec.genes),
        "structure_selected_from": "the module's full response matrix at proposal time",
        "loop": {"max_iters": args.max_iters, "n_starts": args.n_starts,
                 "iterations_recorded": len(result.history)},
        "best": {"spec": json.loads(spec.to_json())},
    }
    with open(args.out, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"{args.module}/{args.condition}: {len(spec.edges)} edges over "
          f"{len(spec.genes)} genes -> {args.out}")
    for h in result.history:
        print(f"  iter loss {h.get('loss')}  structural residuals {h.get('n_structural')}")


if __name__ == "__main__":
    main()
