import numpy as np

from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term
from mmc.loop import anonymize as anon

GENES = ["GATA3", "STAT6", "TBX21", "IL2", "IL2RB", "IL4"]


def _spec():
    return ModelSpec(
        genes=list(GENES),
        edges=[Edge(regulator="STAT6", target="GATA3", sign=1),
               Edge(regulator="TBX21", target="GATA3", sign=-1),
               Edge(regulator="GATA3", target="IL4", sign=1),
               Edge(regulator="IL2RB", target="IL2", sign=-1)],
        rules={"GATA3": Rule(terms=[Term(regulators=["STAT6", "TBX21"],
                                         signs={"TBX21": -1})]),
               "IL4": Rule(terms=[Term(regulators=["GATA3"])]),
               "IL2": Rule(terms=[Term(regulators=["IL2RB"])])},
    )


def test_alias_map_is_a_bijection():
    m = anon.make_alias_map(GENES, seed=0)
    assert set(m.to_alias) == set(GENES)
    assert len(set(m.to_alias.values())) == len(GENES)
    for g, a in m.to_alias.items():
        assert m.to_symbol[a] == g


def test_alias_assignment_is_randomised_across_seeds():
    """A fixed assignment would make the alias itself a name after a few runs."""
    a = anon.make_alias_map(GENES, seed=1).to_alias
    b = anon.make_alias_map(GENES, seed=2).to_alias
    assert a != b
    assert anon.make_alias_map(GENES, seed=1).to_alias == a


def test_substitution_does_not_corrupt_an_overlapping_symbol():
    """The failure the longest-first rule exists to prevent: IL2 inside IL2RB."""
    m = anon.make_alias_map(GENES, seed=0)
    out = anon.anonymise_text("IL2RB represses IL2.", m)
    assert out == f"{m.to_alias['IL2RB']} represses {m.to_alias['IL2']}."
    # the naive failure would leave a dangling "RB" behind the IL2 alias
    assert "RB" not in out
    assert anon.redaction_violations(out, GENES) == []


def test_anonymised_text_carries_no_symbol():
    m = anon.make_alias_map(GENES, seed=3)
    text = "GATA3 activates IL4 and IL5; TBX21 represses GATA3. STAT6 is upstream."
    assert anon.redaction_violations(anon.anonymise_text(text, m), GENES) == []


def test_spec_anonymisation_preserves_structure_and_round_trips():
    m = anon.make_alias_map(GENES, seed=4)
    spec = _spec()
    hidden = anon.anonymise_spec(spec, m)
    assert anon.redaction_violations(hidden.to_json(), GENES) == []
    assert len(hidden.edges) == len(spec.edges)
    back = anon.deanonymise_spec(hidden, m)
    assert {(e.regulator, e.target, e.sign) for e in back.edges} == \
           {(e.regulator, e.target, e.sign) for e in spec.edges}
    assert set(back.rules) == set(spec.rules)
    assert back.rules["GATA3"].terms[0].signs == {"TBX21": -1}


def test_the_audit_catches_a_symbol_that_survived():
    """The audit must fail loudly rather than pass a leaky prompt."""
    assert anon.redaction_violations("the regulator GATA3 acts here", GENES) == ["GATA3"]
    assert anon.redaction_violations("G01 acts on G02", GENES) == []


def test_the_audit_ignores_symbols_embedded_in_longer_words():
    assert anon.redaction_violations("XGATA3Y and GATA3B", GENES) == []


def test_full_prompt_surface_is_clean_after_anonymisation():
    """Covers all four channels the model is shown: gene list, context, spec, residuals."""
    m = anon.make_alias_map(GENES, seed=5)
    surface = anon.assemble_prompt_surface(
        genes=[m.to_alias[g] for g in GENES],
        context=anon.anonymise_text("GATA3 drives IL4 in CD4+ T cells; TBX21 opposes it.", m),
        spec=anon.anonymise_spec(_spec(), m),
        residual_summary=anon.anonymise_text(
            "wrong_sign: knockdown of STAT6 lowers IL2RB but the model raises it", m),
    )
    assert anon.redaction_violations(surface, GENES) == []


def test_an_unanonymised_channel_is_detected():
    """If any single channel is missed the audit must still fail; this is the
    regression guard for the arm being defeated silently at the repair step."""
    m = anon.make_alias_map(GENES, seed=6)
    surface = anon.assemble_prompt_surface(
        genes=[m.to_alias[g] for g in GENES],
        context=anon.anonymise_text("GATA3 drives IL4.", m),
        spec=_spec(),                                   # spec left in real symbols
        residual_summary="",
    )
    assert anon.redaction_violations(surface, GENES)


def test_label_permutation_destroys_pairing_but_keeps_marginals():
    obs = np.arange(24, dtype=float).reshape(6, 4)
    de = np.zeros((6, 4), bool)
    de[:, 0] = True
    p_obs, p_de = anon.permute_perturbation_labels(obs, de, seed=0)
    assert p_obs.shape == obs.shape
    assert sorted(p_obs.sum(axis=1).tolist()) == sorted(obs.sum(axis=1).tolist())
    assert not np.array_equal(p_obs, obs)
    assert p_de.shape == de.shape


def test_label_permutation_is_a_derangement():
    obs = np.arange(40, dtype=float).reshape(8, 5)
    de = np.ones((8, 5), bool)
    p_obs, _ = anon.permute_perturbation_labels(obs, de, seed=2)
    for i in range(8):
        assert not np.array_equal(p_obs[i], obs[i]), f"row {i} kept its own response"


def test_label_permutation_handles_a_degenerate_module():
    obs = np.ones((1, 3))
    de = np.ones((1, 3), bool)
    p_obs, p_de = anon.permute_perturbation_labels(obs, de, seed=0)
    assert p_obs.shape == (1, 3) and p_de.shape == (1, 3)


def test_empty_gene_set_is_handled():
    m = anon.make_alias_map([], seed=0)
    assert m.to_alias == {}
    assert anon.anonymise_text("nothing here", m) == "nothing here"
    assert anon.redaction_violations("nothing here", []) == []
