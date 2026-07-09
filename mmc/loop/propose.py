"""Propose or edit structure via mmc.shared.llm.

The reasoning step sees only the module gene list, the biological context, the
current ModelSpec, and the structural residual summary. It never sees held-out or
test-state data (Tier A) or the Tier-B held-out subset. Every proposal and its
rationale is logged to the trace. Implemented in Step 3.
"""
raise NotImplementedError
