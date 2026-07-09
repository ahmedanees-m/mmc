"""The fit-vs-structure diagnostic gate.

A residual is labelled structural only if no seed recovers its sign or drives it
under tolerance; otherwise it is parametric and withheld, since it is the
optimizer's to resolve, not the structure's. Sign errors are structural; magnitude
errors are parametric. Only structural residuals and the convergence statistics
reach the residual reader.

Contract: diagnose(spec, train, n_starts) -> (best_fit,
{(pert, gene): 'structural' | 'parametric'}, convergence_stats). Implemented in
Step 2, where it labels an injected parametric failure as parametric.
"""
raise NotImplementedError
