"""Shared utilities: receipts and signing, the model-evolution trace, module registration."""
import json

from mmc.data import module_extract
from mmc.grammar.model_spec import Edge, ModelSpec, Rule, Term
from mmc.loop.propose import model_hash
from mmc.shared import receipts
from mmc.shared.trace import Step, Trace


def test_receipt_and_gated():
    r = receipts.Receipt(value=1.5, source="fit", computation="c", query="q")
    g = receipts.Gated(value=1.5, receipt=r)
    assert g.receipt.source == "fit"
    assert r.sig is None


def test_sign_is_deterministic_and_sensitive():
    a = receipts.sign(b"payload", b"key")
    assert a == receipts.sign(b"payload", b"key")
    assert len(a) == 64
    assert receipts.sign(b"other", b"key") != a


def test_key_returns_bytes():
    assert isinstance(receipts._key(), bytes)


def test_trace_log_and_save(tmp_path):
    t = Trace(module="M")
    t.log(Step(kind="propose", model_hash="abc", rationale="r", edit=None, train_score=0.5))
    path = tmp_path / "trace.json"
    t.save(path)
    d = json.loads(path.read_text())
    assert d["module"] == "M"
    assert d["steps"][0]["kind"] == "propose"


def test_register_and_model_genes():
    module_extract.register_module("TESTMOD", ["A", "B"], ["B", "C"])
    assert module_extract.model_genes("TESTMOD") == ["A", "B", "C"]
    assert module_extract.regulators("TESTMOD") == ["A", "B"]


def test_model_hash():
    spec = ModelSpec(genes=["A", "B"],
                     edges=[Edge(regulator="A", target="B", sign=1)],
                     rules={"B": Rule(terms=[Term(regulators=["A"])])})
    h = model_hash(spec)
    assert isinstance(h, str)
    assert len(h) == 10
