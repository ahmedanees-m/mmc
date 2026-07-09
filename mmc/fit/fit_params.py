"""Multi-start parameter fitting (CMA-ES by default, diffrax optional).

Given the structure, fit the continuous parameters to the training-split knockdown
responses only, from at least 16 starts, with bounded parameters for stable
integration and regularization. The multiple starts let the diagnostic gate
separate a fit failure from a structure failure. The optimizer searches over basal,
decay, and per-term production, weight, and threshold values, minimizing the
predicted-versus-observed knockdown delta on the training split, and returns every
seed for diagnose().

Contract: fit(spec, train_deltas) -> (best_params, all_fits). Implemented in Step 2.
"""
raise NotImplementedError


