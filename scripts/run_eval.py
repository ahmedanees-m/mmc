"""Run the Step 5 in-context gate and the Step 6 two-tier transfer test on a frozen
model, then score the conserved and rewired decomposition against the Stage-0 scaffold.

The frozen model is loaded from disk (see freeze_model.py); the Tier B adaptation runs
the loop on the discovery subset only, starting from the frozen structure, so it makes
model calls and needs ANTHROPIC_API_KEY and MMC_ZHU_STORE. Results are written to JSON.
"""
from __future__ import annotations

import json
import sys

from mmc.data.splits import build
from mmc.eval import conserved_rewired, evaluate
from mmc.loop.propose import model_hash
from mmc.loop.run import discover


def _fmt(x) -> str:
    if x is None or (isinstance(x, float) and x != x):
        return "n/a"
    return f"{x:.3f}"


def _table(scores: dict) -> None:
    print(f"    {'method':12s} {'Pearson(all)':>14s} {'sign-acc(DE)':>14s}")
    for k, sc in scores.items():
        pooled = sc.get("pooled") or {}
        print(f"    {k:12s} {_fmt(pooled.get('pearson')):>14s} "
              f"{_fmt(pooled.get('sign_accuracy')):>14s}")


def main() -> None:
    module = sys.argv[1] if len(sys.argv) > 1 else "TCR_signalosome"
    frozen_path = sys.argv[2] if len(sys.argv) > 2 else "/app/frozen_model.json"
    out_path = sys.argv[3] if len(sys.argv) > 3 else "/app/eval_results.json"

    s = build(module)
    frozen_spec, _ = evaluate.load_frozen(frozen_path)

    print("=== in-sample check (frozen model on its own training perturbations) ===")
    ins = evaluate.insample(module, frozen_path)
    _table(ins["scores"])

    print("\n=== Step 5: in-context gate (train state, held-out perturbations) ===")
    gate = evaluate.incontext_gate(module, frozen_path)
    _table(gate["scores"])
    print(f"    gate passed: {gate['passed_gate']}")

    print("\n=== Step 6 Tier A: strict transfer (frozen model to every test-state perturbation) ===")
    a = evaluate.tier_a(module, frozen_path)
    _table(a["scores"])
    print(f"    beats mean: {a['beats_mean']}, beats persistence: {a['beats_persistence']}")

    print("\n=== Step 6 Tier B: adapt on the discovery subset, score the held-out subset ===")
    adapt = discover(module, s.test_state, perts=list(s.tierB_discovery),
                     start_spec=frozen_spec, max_iters=2, n_starts=8, max_iter=50)
    adapted = adapt.ensemble.best()
    print(f"    adapted {model_hash(adapted.spec)} from frozen {model_hash(frozen_spec)}")
    b = evaluate.tier_b(module, adapted.spec, adapted.params)
    _table(b["scores"])
    print(f"    beats consensus: {b['beats_consensus']}, beats mean: {b['beats_mean']}")

    print("\n=== conserved and rewired decomposition against the Stage-0 scaffold ===")
    cr = conserved_rewired.score(frozen_spec, adapted.spec, module)
    print(f"    evaluated edges {cr['evaluated_edges']} (coverage {_fmt(cr['coverage'])})")
    print(f"    rewiring precision {_fmt(cr['precision'])} "
          f"recall {_fmt(cr['recall'])} f1 {_fmt(cr['f1'])}")
    print(f"    predicted rewired: {cr['predicted_rewired']}")
    print(f"    true rewired (evaluated): {cr['true_rewired']}")

    with open(out_path, "w") as f:
        json.dump({"gate": gate, "tier_a": a, "tier_b": b, "conserved_rewired": cr},
                  f, indent=2, default=str)
    print(f"\nsaved {out_path}")


if __name__ == "__main__":
    main()
