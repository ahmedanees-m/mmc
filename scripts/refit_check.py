"""Re-fit a frozen structure's parameters at a large optimizer budget and report the
in-sample fit, to separate an optimizer-budget shortfall from a structure-capacity limit.

If the in-sample Pearson and sign accuracy rise sharply with more iterations, the earlier
poor fit was under-optimization on a structure that can represent the data. If they stay
flat, the structure (or the grammar) cannot fit this module, which is a deeper finding.
"""
from __future__ import annotations

import sys

from mmc.data import module_extract
from mmc.data.splits import build
from mmc.eval import metrics
from mmc.eval.evaluate import load_frozen, model_predict
from mmc.fit.fit_params import multi_fit


def _fit_at(spec, observed, genes, train_perts, n_starts, max_iter) -> None:
    fits = multi_fit(spec, observed, n_starts=n_starts, max_iter=max_iter)
    best = fits[0]
    pred = model_predict(spec, best["params"], train_perts, genes)
    pooled = metrics.score_set(pred, observed, genes, 0.5)["pooled"]
    r = pooled["pearson"]
    s = pooled["sign_accuracy"]
    print(f"  max_iter {max_iter:>4d} n_starts {n_starts:>2d}: training loss {best['loss']:.4f}  "
          f"in-sample Pearson {r if r is None else round(r, 3)}  "
          f"sign-acc {s if s is None else round(s, 3)}")


def main() -> None:
    module = sys.argv[1] if len(sys.argv) > 1 else "CD4_TF_network"
    path = sys.argv[2] if len(sys.argv) > 2 else "/app/frozen_TFnet.json"
    spec, _ = load_frozen(path)
    s = build(module)
    genes = module_extract.model_genes(module)
    train = list(s.train_perts)
    observed = {p: d for p, d in module_extract.observed_deltas(module, s.train_state).items()
                if p in set(train)}

    print(f"refit isolation on {module}: {len(spec.genes)} genes, {len(spec.edges)} edges, "
          f"{len(train)} training perturbations")
    for max_iter in (110, 300, 600):
        _fit_at(spec, observed, genes, train, 12, max_iter)


if __name__ == "__main__":
    main()
