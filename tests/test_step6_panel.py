import importlib.util
from pathlib import Path

import pytest

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "step6_panel.py"
_spec = importlib.util.spec_from_file_location("step6_panel", _PATH)
step6 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step6)


def test_wilson_reproduces_the_published_corpus_intervals():
    """The record reports 0 of 76 novel hypotheses validated, Wilson 95% [0, 4.8%],
    and 0 of 9 module-conditions, [0, 30%]. The implementation must reproduce both or
    the panel's pooled interval is not comparable to the existing figure."""
    lo, hi = step6.wilson(0, 76)
    assert lo == pytest.approx(0.0, abs=1e-9)
    assert hi == pytest.approx(0.048, abs=0.001)

    lo, hi = step6.wilson(0, 9)
    assert lo == pytest.approx(0.0, abs=1e-9)
    assert hi == pytest.approx(0.299, abs=0.002)


def test_wilson_is_symmetric_for_a_half_rate():
    lo, hi = step6.wilson(5, 10)
    assert lo + hi == pytest.approx(1.0, abs=1e-9)


def test_wilson_bounds_stay_inside_the_unit_interval():
    for k, n in ((0, 1), (1, 1), (0, 3), (3, 3), (1, 200)):
        lo, hi = step6.wilson(k, n)
        assert 0.0 <= lo <= hi <= 1.0


def test_wilson_of_an_empty_sample_is_undefined_not_a_crash():
    lo, hi = step6.wilson(0, 0)
    assert lo != lo and hi != hi          # NaN


def test_wilson_narrows_as_the_sample_grows():
    _, hi_small = step6.wilson(0, 10)
    _, hi_large = step6.wilson(0, 1000)
    assert hi_large < hi_small


def test_the_panel_matches_the_locked_amendment():
    """The identifiers are fixed in PREREG_v4 amendment A6; a silent edit here would
    make the panel a different experiment from the pre-registered one."""
    assert step6.PANEL == [
        ("claude-opus-4-8", "anthropic"),
        ("claude-sonnet-5", "anthropic"),
        ("openai/gpt-oss-120b", "ngc"),
        ("meta/llama-3.1-70b-instruct", "ngc"),
        ("nvidia/nemotron-3-ultra-550b-a55b", "ngc"),
        ("z-ai/glm-5.2", "ngc"),
    ]
    families = {m.split("/")[0] for m, p in step6.PANEL if p == "ngc"}
    assert len(families) == 4, "the NGC slots must span four distinct families"
