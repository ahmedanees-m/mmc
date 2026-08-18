"""Diagnostics for the section 8.1 held-out modules, computed before their ceilings exist.

The prospective test predicts five ceilings from diagnostics alone and commits the
predictions before the modules are run. The diagnostics are properties of the observed
response matrix and section 3 defines them as computed before any modelling, so producing
them for a held-out module reveals nothing about its ceiling: the held-out quantity is the
comparator result, which this does not touch.

Kept separate from the comparator so it can run while those modules are still blocked.

    python scripts/step7_holdout_diagnostics.py --holdout prereg/step7_holdout.json \
        --defs /work/step7_defs --condition Stim8hr --out results/step7_holdout_diag.json
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from mmc.data import module_data, module_extract
from mmc.eval import identifiability


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--holdout", required=True)
    ap.add_argument("--defs", required=True, help="directory of module definition JSON")
    ap.add_argument("--condition", default="Stim8hr")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.holdout) as f:
        names = json.load(f)["holdout"]

    rows = {}
    for name in names:
        with open(os.path.join(args.defs, name + ".json")) as f:
            d = json.load(f)
        module_extract.register_module(name, d["regulators"], d["targets"])
        mod = module_data.build_module_data(name, args.condition)
        summary = module_data.module_summary(mod)
        diag = identifiability.diagnostics(mod)
        rows[name] = {"module_summary": summary,
                      "diagnostics": {k: v for k, v in diag.items() if k != "effect_sizes"},
                      "effect_sizes": diag.get("effect_sizes")}
        print(f"{name:<24} genes={summary['n_genes']:<3} perts={summary['n_perts']:<3} "
              f"folds={summary['n_folds_with_de']:<3} "
              f"erank={diag['effective_rank']:.2f} "
              f"pc1={diag['leading_pc_fraction']:.3f} "
              f"spec={diag['perturbation_specific_ratio']:.3f}")

    p = Path(args.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as f:
        json.dump({"condition": args.condition, "modules": rows}, f, indent=2, default=float)
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
