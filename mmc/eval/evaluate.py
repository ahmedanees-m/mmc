"""Two-tier evaluation driver.

Tier A (strict transfer): train-state-frozen ensemble + WT test-state ONLY ->
predict test-state; vs mean/linear/State.
Tier B (few-shot rewiring discovery): + pre-registered discovery subset ->
predict disjoint held-out subset; vs linear/State/consensus.
Freeze before reading. Enforce leakage boundaries.  DoD: Step 6.
"""
raise NotImplementedError
