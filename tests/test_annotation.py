import textwrap

import numpy as np
import pytest

from mmc.eval import annotation as ann

HEADER = "source\ttarget\tsource_genesymbol\ttarget_genesymbol\tis_stimulation\tis_inhibition\n"
ROWS = [
    ("P1", "P2", "GATA3", "IL5", "True", "False"),
    ("P3", "P4", "TBX21", "IL5", "False", "True"),
    ("P5", "P6", "STAT6", "GATA3", "True", "False"),
    ("P7", "P8", "FOO", "BAR", "True", "False"),
    ("P9", "PA", "AMBIG", "IL5", "True", "True"),
    ("PB", "PC", "SELF", "SELF", "True", "False"),
]


@pytest.fixture
def regulon(tmp_path):
    p = tmp_path / "collectri.tsv"
    p.write_text(HEADER + "".join("\t".join(r) + "\n" for r in ROWS), encoding="utf-8")
    return ann.load_collectri(str(p))


def test_loader_reads_symbols_and_signs(regulon):
    assert ("GATA3", "IL5") in regulon.edges
    assert regulon.signed[("GATA3", "IL5")] == 1
    assert regulon.signed[("TBX21", "IL5")] == -1


def test_loader_drops_self_edges(regulon):
    assert ("SELF", "SELF") not in regulon.edges


def test_an_ambiguous_sign_is_unsigned_but_still_an_edge(regulon):
    """Flagged both stimulating and inhibiting carries no usable sign."""
    assert ("AMBIG", "IL5") in regulon.edges
    assert ("AMBIG", "IL5") not in regulon.signed


def test_restriction_keeps_only_edges_inside_the_module(regulon):
    genes = ["GATA3", "STAT6", "IL5"]
    r = regulon.restrict(genes)
    assert r.edges == {("GATA3", "IL5"), ("STAT6", "GATA3")}
    assert ("FOO", "BAR") not in r.edges


def test_coverage_reports_what_the_annotation_speaks_to(regulon):
    cov = regulon.coverage(["GATA3", "STAT6", "IL5", "ORPHAN"])
    assert cov["n_module_genes"] == 4
    assert cov["n_annotated_edges_within_module"] == 2
    assert cov["fraction_genes_covered"] == pytest.approx(0.75)


def test_scoring_a_perfect_unsigned_edge_set(regulon):
    genes = ["GATA3", "STAT6", "IL5"]
    edges = [("GATA3", "IL5", 1), ("STAT6", "GATA3", 1)]
    s = ann.score_edges(edges, regulon, genes)
    assert s["precision"] == pytest.approx(1.0)
    assert s["recall"] == pytest.approx(1.0)
    assert s["jaccard"] == pytest.approx(1.0)


def test_scoring_penalises_a_wrong_sign_only_in_signed_mode(regulon):
    genes = ["GATA3", "TBX21", "IL5"]
    edges = [("TBX21", "IL5", 1)]                 # annotation says this is inhibitory
    unsigned = ann.score_edges(edges, regulon, genes, use_sign=False)
    signed = ann.score_edges(edges, regulon, genes, use_sign=True)
    assert unsigned["precision"] == pytest.approx(1.0)
    assert signed["precision"] == pytest.approx(0.0)


def test_recall_is_against_what_was_recoverable_not_the_whole_regulon(regulon):
    """A module cannot recover edges between genes it does not contain."""
    genes = ["GATA3", "IL5"]
    s = ann.score_edges([("GATA3", "IL5", 1)], regulon, genes)
    assert s["n_annotated_within_module"] == 1
    assert s["recall"] == pytest.approx(1.0)


def test_an_empty_annotation_gives_undefined_recall_not_zero(regulon):
    genes = ["ORPHAN1", "ORPHAN2"]
    s = ann.score_edges([("ORPHAN1", "ORPHAN2", 1)], regulon, genes)
    assert s["n_annotated_within_module"] == 0
    assert np.isnan(s["recall"])


def test_enrichment_separates_a_real_hit_from_chance(regulon):
    genes = ["GATA3", "STAT6", "IL5", "A", "B", "C", "D", "E"]
    good = ann.enrichment([("GATA3", "IL5", 1), ("STAT6", "GATA3", 1)],
                          regulon, genes, n_draws=500)
    bad = ann.enrichment([("A", "B", 1), ("C", "D", 1)], regulon, genes, n_draws=500)
    assert good["observed_hits"] == 2
    assert bad["observed_hits"] == 0
    assert good["p"] < bad["p"]
    assert good["ratio"] > 1.0


def test_enrichment_handles_an_empty_edge_set(regulon):
    out = ann.enrichment([], regulon, ["GATA3", "IL5"], n_draws=50)
    assert np.isnan(out["chance_precision"])
