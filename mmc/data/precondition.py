"""Stage-0 precondition test (regulator + edge level).

For each candidate module x condition pair: pull per-regulator->target effects
(shared.store.module_effects); classify each edge conserved / rewired /
untestable; compute conservation & rewiring fractions. Pass = conservation>=0.5,
rewiring>0, >=N_min testable edges. Produces the Step-0 scaffold + module pick.
DoD: Step 0. (A regulator-breadth fragment already ran on the suppl table;
this is the edge-level version against the full DE store.)
"""
raise NotImplementedError
