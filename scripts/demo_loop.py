"""Run the discovery loop on a module and print the evolution trace, the training-fit
history, and the frozen ensemble with per-edge agreement. Requires MMC_ZHU_STORE and
ANTHROPIC_API_KEY.
"""
from __future__ import annotations

import sys

from mmc.loop.run import discover


def main() -> None:
    module = sys.argv[1] if len(sys.argv) > 1 else "TCR_core"
    condition = sys.argv[2] if len(sys.argv) > 2 else "Stim8hr"
    result = discover(module, condition, max_iters=2, n_starts=4, max_iter=30)

    print(f"=== discovery loop: {module} in {condition} ===\n")
    print("training-fit history (loss should improve or plateau):")
    for h in result.history:
        print(f"  {h['hash']}  loss {h['loss']:.4f}  structural residuals {h['n_structural']}")

    print("\nmodel-evolution trace:")
    for s in result.trace.steps:
        score = f" loss {s.train_score:.4f}" if s.train_score is not None else ""
        edit = f"  edit={s.edit}" if s.edit else ""
        print(f"  [{s.kind}] {s.model_hash}{score}  {s.rationale}{edit}")

    best = result.ensemble.best()
    print(f"\nfrozen ensemble: {len(result.ensemble.members)} structure(s), "
          f"best training loss {best.loss:.4f}")
    print("edge agreement across the ensemble:")
    for edge, frac in sorted(result.ensemble.edge_agreement().items(), key=lambda kv: -kv[1]):
        regulator, target, sign = edge
        print(f"  {regulator} -> {target} ({sign:+d}): {frac:.2f}")


if __name__ == "__main__":
    main()
