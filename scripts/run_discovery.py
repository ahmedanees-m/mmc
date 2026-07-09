"""v3 Part 3: run the discovery loop as an instrument on an immune module.

With the structural backend the loop proposes a structure, fits it in seconds, reads the
genuine structural residuals through the fixed gate, reasons about them, and revises the
structure. This prints the training-fit history, the engineer's proposal-and-rationale
trace (including the falsified steps), and the frozen structure, so the discovery
candidates and the reasoning behaviour can be read off. Requires the model client.
"""
from __future__ import annotations

import sys

from mmc.compile import structural
from mmc.data import module_extract
from mmc.eval import metrics
from mmc.loop.propose import model_hash
from mmc.loop.run import discover


def _insample(spec, params, module, condition) -> dict:
    genes = module_extract.model_genes(module)
    observed = module_extract.observed_deltas(module, condition)
    gi = {g: i for i, g in enumerate(spec.genes)}
    wt = structural.steady_state(spec, params)
    pred = {}
    for p in observed:
        if p in gi:
            d = structural.knockdown(spec, params, p, wt=wt)
            pred[p] = {g: float(d[gi[g]]) for g in genes if g in gi and g != p}
    return metrics.score_set(pred, observed, genes, 0.5)["pooled"]


def main() -> None:
    module = sys.argv[1] if len(sys.argv) > 1 else "CD4_lineage_TFs"
    condition = sys.argv[2] if len(sys.argv) > 2 else "Stim8hr"
    max_iters = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    result = discover(module, condition, backend="structural",
                      max_iters=max_iters, n_starts=8, max_iter=300)

    print(f"=== discovery loop: {module} in {condition} (structural backend) ===\n")
    print("training-fit history (loss should improve or plateau):")
    for h in result.history:
        print(f"  {h['hash']}  loss {h['loss']:.4f}  structural residuals {h['n_structural']}")

    print("\nengineer trace (proposals, rationales, falsified steps):")
    for s in result.trace.steps:
        score = f" loss {s.train_score:.4f}" if s.train_score is not None else ""
        edit = f"  edit={s.edit}" if s.edit else ""
        print(f"  [{s.kind}] {s.model_hash}{score}  {s.rationale}{edit}")

    best = result.ensemble.best()
    pooled = _insample(best.spec, best.params, module, condition)
    print(f"\nfrozen structure {model_hash(best.spec)}: {len(best.spec.edges)} edges, "
          f"training loss {best.loss:.4f}")
    print(f"in-sample fit: Pearson {round(pooled['pearson'], 3)}  "
          f"sign-acc {round(pooled['sign_accuracy'], 3)}")
    for e in sorted(best.spec.edges, key=lambda e: (e.target, e.regulator)):
        print(f"  {e.regulator} -> {e.target} ({e.sign:+d})")


if __name__ == "__main__":
    main()
