import pytest

from mmc.data import textbook
from mmc.data.precondition import MODULES
from mmc.grammar.model_spec import MAX_REGS_PER_TERM


def _genes(module):
    m = MODULES[module]
    seen, out = set(), []
    for g in list(m["regulators"]) + list(m["targets"]):
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


@pytest.mark.parametrize("module", sorted(textbook.TEXTBOOK_EDGES))
def test_every_recorded_edge_is_well_formed(module):
    for reg, tgt, sign in textbook.TEXTBOOK_EDGES[module]:
        assert reg != tgt, f"{module}: self edge on {reg}"
        assert sign in (1, -1), f"{module}: bad sign on {reg}->{tgt}"


@pytest.mark.parametrize("module", sorted(textbook.TEXTBOOK_EDGES))
def test_no_duplicate_or_contradictory_edges(module):
    pairs = [(r, t) for r, t, _ in textbook.TEXTBOOK_EDGES[module]]
    assert len(pairs) == len(set(pairs)), f"{module}: an edge is listed twice"


@pytest.mark.parametrize("module", sorted(textbook.TEXTBOOK_EDGES))
def test_recorded_structures_respect_the_in_degree_cap(module):
    """The grammar allows at most three regulators per target, so a textbook set
    that names more has to be trimmed when it is written, not when it is compiled."""
    in_deg: dict[str, int] = {}
    for _r, t, _s in textbook.TEXTBOOK_EDGES[module]:
        in_deg[t] = in_deg.get(t, 0) + 1
    over = {t: n for t, n in in_deg.items() if n > MAX_REGS_PER_TERM}
    assert not over, f"{module}: targets over the in-degree cap: {over}"


@pytest.mark.parametrize("module", ["Th2_GATA3", "TCR_signalosome", "CD4_lineage_TFs"])
def test_spec_builds_for_the_registered_modules(module):
    spec = textbook.textbook_spec(module, _genes(module))
    assert spec.edges
    for e in spec.edges:
        assert e.regulator in spec.genes and e.target in spec.genes


def test_th2_carries_the_mutual_antagonism():
    spec = textbook.textbook_spec("Th2_GATA3", _genes("Th2_GATA3"))
    assert spec.edge_sign("GATA3", "TBX21") == -1
    assert spec.edge_sign("TBX21", "GATA3") == -1
    assert spec.edge_sign("GATA3", "IL5") == 1


def test_tcr_cascade_is_ordered_and_terminates_on_il2():
    spec = textbook.textbook_spec("TCR_signalosome", _genes("TCR_signalosome"))
    for a, b in (("CD3E", "ZAP70"), ("ZAP70", "LAT"), ("LAT", "LCP2"),
                 ("PLCG1", "PRKCQ")):
        assert spec.edge_sign(a, b) == 1, f"missing {a}->{b}"
    assert spec.edge_sign("NFKB1", "IL2") == 1


def test_edges_outside_the_gene_set_are_dropped_not_renamed():
    spec = textbook.textbook_spec("Th2_GATA3", ["GATA3", "STAT6", "IL5"])
    assert {(e.regulator, e.target) for e in spec.edges} == {("STAT6", "GATA3"),
                                                             ("GATA3", "IL5")}


def test_coverage_reports_what_the_module_could_not_hold():
    cov = textbook.coverage("Th2_GATA3", ["GATA3", "STAT6", "IL5"])
    assert cov["n_textbook_edges"] == 9
    assert cov["n_edges_in_module"] == 2
    assert "GATA3->IL13" in cov["dropped"]


def test_unknown_module_raises():
    with pytest.raises(KeyError):
        textbook.textbook_spec("NoSuchModule", ["A"])
