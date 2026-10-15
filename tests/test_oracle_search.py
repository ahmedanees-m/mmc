import numpy as np

from mmc.eval import oracle_search as osr
from mmc.eval.holdout import ModuleData
from mmc.grammar.model_spec import MAX_REGS_PER_TERM

GENES = ["A", "B", "C", "D"]
PERTS = ["A", "B", "C"]
TARGET_EDGES = {("A", "B", 1), ("C", "D", -1)}


def _mod():
    rng = np.random.default_rng(0)
    obs = rng.normal(size=(3, 4))
    de = np.zeros((3, 4), bool)
    de[:, 1] = True
    return ModuleData(GENES, PERTS, obs, de, None)


def _scorer(calls=None):
    """Score by overlap with a planted edge set, so the search has a known optimum."""
    def score_many(edge_sets):
        out = []
        for es in edge_sets:
            if calls is not None:
                calls.append(frozenset(es))
            hit = len(TARGET_EDGES & set(es))
            penalty = 0.01 * max(0, len(es) - len(TARGET_EDGES))
            score = hit / len(TARGET_EDGES) - penalty
            out.append((score, 1.0 - score))
        return out
    return score_many


def test_spec_from_edges_builds_a_valid_spec():
    spec = osr.spec_from_edges(GENES, TARGET_EDGES)
    assert {(e.regulator, e.target, e.sign) for e in spec.edges} == TARGET_EDGES
    for tgt, rule in spec.rules.items():
        for term in rule.terms:
            for r in term.regulators:
                assert any(e.regulator == r and e.target == tgt for e in spec.edges)


def test_spec_from_edges_enforces_the_in_degree_cap():
    edges = {(r, "D", 1) for r in ["A", "B", "C"]} | {("A", "B", 1)}
    spec = osr.spec_from_edges(["A", "B", "C", "D"], edges)
    in_deg: dict[str, int] = {}
    for e in spec.edges:
        in_deg[e.target] = in_deg.get(e.target, 0) + 1
    assert max(in_deg.values()) <= MAX_REGS_PER_TERM


def test_candidate_pool_is_ranked_and_bounded():
    mod = _mod()
    pool = osr.candidate_pool(mod, size=5)
    assert len(pool) == 5
    assert all(r != t for r, t in pool)
    assert all(r in PERTS and t in GENES for r, t in pool)
    # DE entries outrank non-DE ones, so column B should appear early
    assert any(t == "B" for r, t in pool[:3])


def test_all_pairs_excludes_self_edges():
    pairs = osr.all_pairs(_mod())
    assert len(pairs) == len(PERTS) * len(GENES) - len(PERTS)
    assert all(r != t for r, t in pairs)


def test_greedy_recovers_the_planted_structure():
    mod = _mod()
    found = osr.greedy_forward_backward(_scorer(), osr.all_pairs(mod), max_edges=6)
    assert TARGET_EDGES <= found


def test_greedy_stops_adding_once_nothing_helps():
    mod = _mod()
    found = osr.greedy_forward_backward(_scorer(), osr.all_pairs(mod), max_edges=20)
    # the size penalty means the optimum is exactly the planted set
    assert found == TARGET_EDGES


def test_greedy_backward_removes_a_planted_useless_edge():
    mod = _mod()
    junk = ("B", "C", 1)
    found = osr.greedy_forward_backward(_scorer(), osr.all_pairs(mod),
                                        max_edges=6, init={junk})
    assert junk not in found


def test_greedy_records_a_trace():
    mod = _mod()
    trace = osr.SearchTrace()
    osr.greedy_forward_backward(_scorer(), osr.all_pairs(mod), max_edges=4, trace=trace)
    assert len(trace.holdout) == len(trace.edge_sets) == len(trace.train_loss)
    assert len(trace.holdout) > 1
    assert trace.best_index() >= 0
    d = trace.to_dict()
    assert d["n_evaluated"] == len(trace.holdout)


def test_annealing_improves_on_a_deliberately_bad_start():
    mod = _mod()
    start = {("B", "C", 1)}
    best = osr.simulated_annealing(_scorer(), start, osr.all_pairs(mod),
                                   n_steps=400, batch=16, seed=0)
    assert len(TARGET_EDGES & best) >= 1


def test_annealing_can_reach_edges_outside_a_restricted_pool():
    """The property that keeps the greedy pool from silently capping the ceiling."""
    mod = _mod()
    calls: list = []
    restricted = [("A", "B")]
    start = osr.greedy_forward_backward(_scorer(), restricted, max_edges=4)
    osr.simulated_annealing(_scorer(calls), start, osr.all_pairs(mod),
                            n_steps=600, batch=16, seed=1)
    reached = {e for es in calls for e in es}
    assert any(e[:2] == ("C", "D") for e in reached)


def test_annealing_respects_the_edge_cap():
    mod = _mod()
    best = osr.simulated_annealing(_scorer(), set(), osr.all_pairs(mod),
                                   n_steps=400, batch=8, max_edges=3, seed=2)
    assert len(best) <= 3


def test_screening_shortlists_but_does_not_change_what_a_candidate_is_worth():
    """Screening may reorder which candidates are considered; the accepted score is
    always the real objective, so the planted optimum still has to win on it."""
    mod = _mod()
    real_calls: list = []
    screen_calls: list = []

    def screen_many(sets):
        for es in sets:
            screen_calls.append(es)
        # a deliberately noisy proxy that still ranks the planted edges highly
        return [(len(TARGET_EDGES & set(es)) * 0.5, 1.0) for es in sets]

    found = osr.greedy_forward_backward(_scorer(real_calls), osr.all_pairs(mod),
                                        max_edges=6, screen_many=screen_many,
                                        screen_keep=4)
    assert TARGET_EDGES <= found
    # the expensive objective saw far fewer candidates than the screen did
    assert len(real_calls) < len(screen_calls)


def test_screening_is_bypassed_when_the_batch_is_already_small():
    screen_calls: list = []

    def screen_many(sets):
        screen_calls.extend(sets)
        return [(0.0, 1.0) for _ in sets]

    osr.greedy_forward_backward(_scorer(), [("A", "B")], max_edges=2,
                                screen_many=screen_many, screen_keep=50)
    assert screen_calls == []


def test_greedy_starts_from_an_unscoreable_empty_structure():
    """The empty structure scores NaN on a sparse module, because no fold has a DE
    gene to overlap with. Comparing a candidate against NaN never succeeds, so a naive
    loop returns the empty set having accepted nothing."""
    def score_many(edge_sets):
        out = []
        for es in edge_sets:
            if not es:
                out.append((float("nan"), float("inf")))
            else:
                out.append((len(TARGET_EDGES & set(es)) / len(TARGET_EDGES), 1.0))
        return out

    found = osr.greedy_forward_backward(score_many, osr.all_pairs(_mod()), max_edges=4)
    assert found, "greedy returned the empty set from a NaN start"
    assert TARGET_EDGES <= found


def test_greedy_stops_cleanly_when_every_candidate_is_undefined():
    def score_many(edge_sets):
        return [(float("nan"), float("inf")) for _ in edge_sets]

    found = osr.greedy_forward_backward(score_many, osr.all_pairs(_mod()), max_edges=4)
    assert found == set()


def test_annealing_survives_an_all_nan_objective():
    def score_many(edge_sets):
        return [(float("nan"), float("inf")) for _ in edge_sets]

    best = osr.simulated_annealing(score_many, set(), osr.all_pairs(_mod()),
                                   n_steps=60, batch=6, seed=0)
    assert best == set()


def test_greedy_respects_the_round_cap():
    """A scorer that always claims an improvement would loop forever without the cap."""
    calls = {"n": 0}

    def score_many(edge_sets):
        calls["n"] += 1
        return [(float(calls["n"]), 1.0) for _ in edge_sets]

    osr.greedy_forward_backward(score_many, osr.all_pairs(_mod()),
                                max_edges=2, max_rounds=2)
    assert calls["n"] < 100


def test_edges_of_round_trips_through_spec_from_edges():
    spec = osr.spec_from_edges(GENES, TARGET_EDGES)
    assert osr.edges_of(spec) == TARGET_EDGES
