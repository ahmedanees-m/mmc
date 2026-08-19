"""The Step 3 topology families must actually differ from the uniform sampler.

Amendment A14 adds hub, scale-free and modular generators so the calibration sweep varies
how edges are arranged and not only how many there are. A sampler that names a topology but
draws something close to uniform would make that arm a second uniform arm under another
name, and the calibration would look like it had addressed the objection without doing so.
The scale-free sampler was in exactly that state at its first default, so these thresholds
are here to keep any future parameter change honest.
"""
import collections

import numpy as np

from mmc.eval import random_null as rn

GENES = ["G%02d" % i for i in range(28)]
N_EDGES = 48


def _out_degree_share(spec, top=3):
    counts = collections.Counter(e.regulator for e in spec.edges)
    ordered = sorted(counts.values(), reverse=True)
    return sum(ordered[:top]) / sum(ordered)


def _within_block_fraction(spec, seed, n_blocks=4):
    rng = np.random.default_rng(seed)
    order = list(rng.permutation(len(GENES)))
    block = {GENES[gi]: rank % n_blocks for rank, gi in enumerate(order)}
    hits = sum(1 for e in spec.edges if block[e.regulator] == block[e.target])
    return hits / len(spec.edges)


def _mean(fn, stat, seeds=5):
    return sum(stat(fn(np.random.default_rng(s)), s) for s in range(seeds)) / seeds


def test_every_family_draws_the_requested_edge_count():
    for name, fn in rn.TOPOLOGIES.items():
        spec = fn(GENES, GENES, N_EDGES, np.random.default_rng(0))
        assert len(spec.edges) == N_EDGES, name
        assert len({(e.regulator, e.target) for e in spec.edges}) == N_EDGES, name
        assert all(e.regulator != e.target for e in spec.edges), name


def test_hub_and_scale_free_concentrate_out_degree():
    uniform = _mean(lambda r: rn.sample_spec(GENES, GENES, N_EDGES, r),
                    lambda s, _: _out_degree_share(s))
    hub = _mean(lambda r: rn.sample_spec_hub(GENES, GENES, N_EDGES, r),
                lambda s, _: _out_degree_share(s))
    free = _mean(lambda r: rn.sample_spec_scale_free(GENES, GENES, N_EDGES, r),
                 lambda s, _: _out_degree_share(s))
    # both must be clearly above uniform, or the arm is not testing a hub topology
    assert hub > uniform + 0.15, (hub, uniform)
    assert free > uniform + 0.15, (free, uniform)


def test_modular_concentrates_edges_inside_blocks():
    uniform = _mean(lambda r: rn.sample_spec(GENES, GENES, N_EDGES, r),
                    _within_block_fraction)
    modular = _mean(lambda r: rn.sample_spec_modular(GENES, GENES, N_EDGES, r),
                    _within_block_fraction)
    assert uniform < 0.40, uniform          # near the 0.25 chance level at four blocks
    assert modular > 0.60, modular


def test_in_degree_cap_is_respected_by_every_family():
    for name, fn in rn.TOPOLOGIES.items():
        spec = fn(GENES, GENES, N_EDGES, np.random.default_rng(1))
        by_target = collections.Counter(e.target for e in spec.edges)
        assert max(by_target.values()) <= rn.MAX_REGS_PER_TERM, name
        for target, rule in spec.rules.items():
            assert len(rule.terms) == 1, name
            assert len(rule.terms[0].regulators) <= rn.MAX_REGS_PER_TERM, name
