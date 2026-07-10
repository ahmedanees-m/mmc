"""Moonen 2026 CRE-gene loader.

Used only for orthogonal edge validation: check high-agreement model edges against
variant->CRE->gene links. The CRE-gene table availability is checked first; the
functionally-significant subset may be absent from the supplementary, in which case
the prioritized candidate list is used.
"""
# TODO: expose corroborates(regulator, target) -> bool | None against the store.
