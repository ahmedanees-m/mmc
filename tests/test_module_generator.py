import numpy as np

from mmc.data import module_generator as mg
from mmc.eval.annotation import Regulon


def _regulon():
    edges = {("TF1", f"T{i}") for i in range(10)}
    edges |= {("TF2", f"T{i}") for i in range(3)}          # too few targets
    edges |= {("TF3", f"U{i}") for i in range(8)}          # targets not measured
    return Regulon(edges, {}, "test")


MEASURED = {"TF1", "TF2"} | {f"T{i}" for i in range(10)}
PERTURBED = {"TF1", "T0", "T1", "T2"}


def test_regulon_candidates_need_enough_measured_targets():
    cands = mg.regulon_candidates(_regulon(), MEASURED, PERTURBED, min_targets=6)
    names = {c.name for c in cands}
    assert "regulon_TF1" in names
    assert "regulon_TF2" not in names, "3 targets is below the minimum"
    assert "regulon_TF3" not in names, "its targets are not measured"


def test_regulon_candidates_need_perturbed_genes():
    """A fold is a held-out perturbation, so a module of unperturbed genes has none."""
    cands = mg.regulon_candidates(_regulon(), MEASURED, set(), min_targets=6)
    assert cands == []


def test_regulon_candidate_respects_the_gene_cap_and_keeps_perturbed_targets():
    big = Regulon({("TF", f"G{i}") for i in range(80)}, {}, "test")
    measured = {"TF"} | {f"G{i}" for i in range(80)}
    perturbed = {f"G{i}" for i in range(70, 80)}
    cands = mg.regulon_candidates(big, measured, perturbed, min_targets=6, max_genes=20)
    assert cands and len(cands[0].genes) <= 20
    kept = set(cands[0].genes)
    assert len(kept & perturbed) >= 5, "perturbed targets must survive the trim"


def test_coresponse_clusters_correlated_genes_together():
    rng = np.random.default_rng(0)
    n_perts = 60
    program = rng.normal(size=n_perts)
    cols = [program * rng.uniform(0.8, 1.2) + 0.01 * rng.normal(size=n_perts)
            for _ in range(6)]                            # six genes on one program
    cols += [rng.normal(size=n_perts) for _ in range(6)]  # six independent genes
    x = np.stack(cols, axis=1)
    genes = [f"P{i}" for i in range(6)] + [f"I{i}" for i in range(6)]

    cands = mg.coresponse_candidates(genes, x, size=6, n_modules=1)
    assert cands
    members = set(cands[0].genes)
    on_program = {g for g in members if g.startswith("P")}
    assert len(on_program) >= 5, f"expected the shared program to cluster, got {members}"


def test_clusters_are_seeded_by_co_variation_not_by_variance():
    """Regression test for a real defect.

    Seeding on a gene's own variance picks whichever gene is noisiest, which is the
    opposite of what a co-response module is. This fixture is built so the independent
    genes have the larger standard deviations while the shared-program genes correlate
    almost perfectly with each other; variance-seeding builds its cluster around an
    independent gene and recovers one program gene out of six.
    """
    rng = np.random.default_rng(0)
    n = 40
    program = rng.normal(size=n) * 0.5                    # deliberately low variance
    cols = [program * rng.uniform(0.9, 1.1) + 0.01 * rng.normal(size=n)
            for _ in range(6)]
    cols += [rng.normal(size=n) * 2.0 for _ in range(6)]  # high variance, independent
    x = np.stack(cols, axis=1)
    genes = [f"P{i}" for i in range(6)] + [f"I{i}" for i in range(6)]

    sd = (x - x.mean(axis=0)).std(axis=0, ddof=1)
    assert sd[6:].min() > sd[:6].max(), "fixture must make the noisy genes highest-variance"

    members = set(mg.coresponse_candidates(genes, x, size=6, n_modules=1)[0].genes)
    on_program = {g for g in members if g.startswith("P")}
    assert len(on_program) >= 5, (
        f"the cluster should follow the shared program, not the variance; got {members}")


def test_coresponse_handles_a_module_too_small_to_cluster():
    assert mg.coresponse_candidates(["A", "B"], np.ones((2, 2)), size=6) == []


def test_screen_records_every_failure_reason():
    c = mg.Candidate(name="m", source="regulon", genes=["A", "B"],
                     n_perts=2, n_de_entries=3, n_scoreable_folds=1)
    mg.screen(c)
    assert c.passed is False
    assert len(c.reasons) == 3
    assert any("perturbations" in r for r in c.reasons)
    assert any("DE entries" in r for r in c.reasons)
    assert any("scoreable folds" in r for r in c.reasons)


def test_screen_passes_a_powered_candidate():
    c = mg.Candidate(name="m", source="regulon", genes=["A"] * 20,
                     n_perts=20, n_de_entries=100, n_scoreable_folds=15)
    mg.screen(c)
    assert c.passed is True and c.reasons == []


def test_the_fold_floor_would_have_caught_th2():
    """Th2_GATA3 ran with 7 perturbations, 5 DE entries and 2 scoreable folds, and
    produced intervals spanning the entire metric. The screen must reject that."""
    c = mg.Candidate(name="Th2_GATA3", source="regulon", genes=["A"] * 7,
                     n_perts=7, n_de_entries=5, n_scoreable_folds=2)
    mg.screen(c)
    assert c.passed is False


def test_characterise_counts_scoreable_folds():
    obs = np.zeros((4, 5))
    de = np.zeros((4, 5), bool)
    de[0, 0] = de[2, 1] = de[2, 3] = True
    out = mg.characterise(["a"] * 5, obs, de)
    assert out["n_perts"] == 4
    assert out["n_de_entries"] == 3
    assert out["n_scoreable_folds"] == 2


def test_summary_reports_failures_not_just_passes():
    cands = [
        mg.screen(mg.Candidate("a", "regulon", ["g"], n_perts=20, n_de_entries=99,
                               n_scoreable_folds=15)),
        mg.screen(mg.Candidate("b", "coresponse", ["g"], n_perts=2, n_de_entries=1,
                               n_scoreable_folds=1)),
    ]
    s = mg.summarise(cands)
    assert s["n_extracted"] == 2 and s["n_passed"] == 1 and s["n_failed"] == 1
    assert s["by_source"] == {"regulon": 1, "coresponse": 1}
    assert s["failure_reasons"]


def test_candidate_serialises():
    c = mg.screen(mg.Candidate("a", "regulon", ["X", "Y"], seed_gene="X",
                               n_perts=9, n_de_entries=40, n_scoreable_folds=9))
    d = c.as_dict()
    assert d["passed"] is True and d["n_genes"] == 2 and d["seed_gene"] == "X"
