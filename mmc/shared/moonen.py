"""Moonen 2026 CRE-gene loader, carried over from the VERDICT Moonen store.

Used only for Step-7 orthogonal edge validation: check MMC's high-agreement
model edges against variant->CRE->gene links. Verify the CRE-gene table
availability first (see the VERDICT integrity note: the functionally-significant
subset may be absent from the supplementary; the prioritized candidate list is
what is available).
"""
# TODO: expose corroborates(regulator, target) -> bool | None against the store.
