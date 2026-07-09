"""Ensemble tracker (fixes identifiability).

Steady-state data is degenerate -- many structures fit equally. Keep the top-k
structures within delta-loss / delta-BIC of the best; track per-edge agreement.
The frozen artifact is the ENSEMBLE + its edge-agreement profile; the
conserved/rewired call is made per-edge across it.  DoD: Step 3.
"""
raise NotImplementedError
