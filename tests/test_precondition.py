"""Precondition edge-classification logic, tested offline with synthetic effects
(no store access). Each tuple is (effect_size, fdr, crossguide_r, n_downstream);
the activity n_downstream is the power gate."""
from mmc.data.precondition import ACTIVE_MIN, _label

HI = ACTIVE_MIN + 500   # clearly active
LO = ACTIVE_MIN - 1     # inactive


def edge(effect, fdr, activity):
    return (effect, fdr, None, activity)


def test_conserved_same_sign_significant_both():
    assert _label(edge(-2.0, 1e-4, HI), edge(-1.5, 1e-3, HI)) == "conserved"


def test_rewired_opposite_sign():
    assert _label(edge(-2.0, 1e-4, HI), edge(+2.0, 1e-4, HI)) == "rewired"


def test_rewired_significant_in_one_state_only():
    assert _label(edge(-2.0, 1e-4, HI), edge(-0.1, 0.9, HI)) == "rewired"


def test_no_effect_when_active_but_never_significant():
    assert _label(edge(0.0, 0.9, HI), edge(0.1, 0.8, HI)) == "no_effect"


def test_untestable_when_a_perturbation_is_inactive():
    assert _label(edge(-2.0, 1e-4, LO), edge(-2.0, 1e-4, HI)) == "untestable"


def test_untestable_when_not_measured():
    assert _label(None, edge(-2.0, 1e-4, HI)) == "untestable"


def test_the_four_curated_modules_are_registered():
    """The paper reports four curated modules; all four must be resolvable by name."""
    from mmc.data.precondition import MODULES
    assert set(MODULES) == {"Th2_GATA3", "TCR_signalosome", "CD4_lineage_TFs",
                            "Cytokine_production"}


def test_cytokine_production_matches_its_deposited_shape():
    """28 genes from 18 regulators and 10 readouts, which is what the results report."""
    from mmc.data import module_extract
    from mmc.data.precondition import MODULES

    m = MODULES["Cytokine_production"]
    assert len(m["regulators"]) == 18
    assert len(m["targets"]) == 28
    assert set(m["regulators"]) < set(m["targets"])
    genes = module_extract.model_genes("Cytokine_production")
    assert len(genes) == 28 and len(set(genes)) == 28
    # the regulators come first, which is the gene order the exported matrices carry
    assert genes[:18] == m["regulators"]


def test_every_curated_module_has_regulators_inside_its_targets():
    from mmc.data import module_extract
    from mmc.data.precondition import MODULES

    for name, m in MODULES.items():
        genes = set(module_extract.model_genes(name))
        assert set(m["regulators"]) <= genes, name
        assert set(m["targets"]) <= genes, name
