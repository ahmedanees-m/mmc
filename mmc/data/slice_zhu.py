"""Slice a module's per-condition knockdown effect table from the store, for
inspection. This is the module effect table the fitting target is assembled from
(see module_extract). Run as a script to print all three conditions.
"""
from __future__ import annotations

from ..shared import store
from . import module_extract


def slice_module(module: str, condition: str):
    """The regulator-to-gene effect table for a module in one condition."""
    genes = module_extract.model_genes(module)
    regs = module_extract.regulators(module)
    return store.module_effects(regs, genes, condition)


def main(module: str = "TCR_signalosome") -> None:
    for condition in store.CONDITIONS:
        df = slice_module(module, condition)
        print(f"\n=== {module}  {condition}  ({len(df)} effects) ===")
        print(df.to_string(index=False))


if __name__ == "__main__":
    import sys

    main(sys.argv[1] if len(sys.argv) > 1 else "TCR_signalosome")
