"""Assemble a module's perturbation-response training data from the store.

For a module and a condition, gather the observed knockdown effect (log2 fold
change) of each candidate regulator on each module gene, at module resolution (the
per-gene shift the store already holds). This is the target the fitter matches:
predicted knockdown delta against observed. The on-target self effect is excluded,
since it is imposed by the perturbation operator rather than predicted.
"""
from __future__ import annotations

from ..shared import store
from .precondition import MODULES

# Modules assembled at runtime register here, so a gene list produced by a screen can be
# run by name without being written into MODULES first. `scripts/assemble_modules.py`
# writes such lists to `data/module_definitions.json`, and `run_comparator.py --module-def`
# registers one. The four curated modules are in MODULES and need neither.
_DYNAMIC: dict[str, dict[str, list[str]]] = {}


def register_module(name: str, regulators: list[str], targets: list[str]) -> None:
    _DYNAMIC[name] = {"regulators": list(regulators), "targets": list(targets)}


def _spec(module: str) -> dict:
    return _DYNAMIC.get(module) or MODULES[module]


def model_genes(module: str) -> list[str]:
    """The genes in the model: the union of regulators and targets, order-stable."""
    m = _spec(module)
    seen: set[str] = set()
    out: list[str] = []
    for g in list(m["regulators"]) + list(m["targets"]):
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def regulators(module: str) -> list[str]:
    return list(_spec(module)["regulators"])


def observed_deltas(module: str, condition: str) -> dict[str, dict[str, float]]:
    """{perturbation: {gene: observed log2 fold change}} for the module in one condition.

    Only module genes are kept, and only knockdowns of module genes (the genes the
    model can perturb). The on-target self effect is dropped.
    """
    genes = model_genes(module)
    regs = [g for g in regulators(module) if g in genes]
    df = store.module_effects(regs, genes, condition)
    out: dict[str, dict[str, float]] = {}
    for _, row in df.iterrows():
        pert, gene = row["perturbation"], row["target_gene"]
        if pert == gene:
            continue
        out.setdefault(pert, {})[gene] = float(row["effect_size"])
    return out
