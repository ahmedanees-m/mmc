"""Residual reader: compact, structured summary of SYSTEMATIC STRUCTURAL failures.

Groups by pattern: wrong-sign / missing-effect / spurious-effect / non-monotone
(the last signals a logic-gate repair). Attaches (perturbation, target, gap) +
fit-confidence. Withholds parametric (magnitude-only) residuals. Token-budgeted.
Contract: summarize(structural_residuals, convergence) -> dict.  DoD: Step 3.
"""
raise NotImplementedError
