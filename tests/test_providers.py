import json

import pytest

from mmc.loop import providers as P


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for var in ("MMC_MODEL", "MMC_PROVIDER", "MMC_BASE_URL",
                "NVIDIA_API_KEY", "NGC_API_KEY", "OPENAI_API_KEY"):
        monkeypatch.delenv(var, raising=False)
    P.STATS.__init__()


def _provider(**kw):
    return P.OpenAICompatProvider("test/model", "https://example.invalid/v1", "k", **kw)


def _reply(content=None, reasoning=None, tokens=11):
    msg = {}
    if content is not None:
        msg["content"] = content
    if reasoning is not None:
        msg["reasoning_content"] = reasoning
    return {"choices": [{"message": msg}], "usage": {"total_tokens": tokens}}


def test_answer_in_content_only():
    assert _provider()._read(_reply(content="hello")) == "hello"


def test_answer_in_reasoning_only_is_not_lost():
    """The bug this guards: reading only `content` scored gpt-oss-120b as empty when
    it had answered, which would have been reported as a model failure."""
    p = _provider()
    assert p._read(_reply(content="", reasoning="the answer")) == "the answer"
    assert P.STATS.empty_content_with_reasoning == 1


def test_both_fields_are_returned_reasoning_first():
    out = _provider()._read(_reply(content='```json\n{"a":1}\n```', reasoning="because"))
    assert out.startswith("because")
    assert '{"a":1}' in out


def test_missing_message_yields_empty_string_not_a_crash():
    assert _provider()._read({"choices": [{}]}) == ""
    assert _provider()._read({}) == ""


def test_successful_call_records_stats(monkeypatch):
    p = _provider()
    monkeypatch.setattr(p, "_post", lambda payload: _reply(content="ok", tokens=42))
    assert p.complete("sys", "user") == "ok"
    assert P.STATS.calls == 1
    assert P.STATS.total_tokens == 42
    assert P.STATS.by_model["test/model"] == 1


def test_the_prompt_is_sent_as_system_plus_user_with_temperature_zero(monkeypatch):
    seen = {}
    p = _provider()

    def fake(payload):
        seen.update(payload)
        return _reply(content="ok")

    monkeypatch.setattr(p, "_post", fake)
    p.complete("SYSTEM", "USER", max_tokens=123)
    assert seen["messages"] == [{"role": "system", "content": "SYSTEM"},
                                {"role": "user", "content": "USER"}]
    assert seen["temperature"] == 0
    assert seen["max_tokens"] == 123
    assert seen["model"] == "test/model"


def test_retries_a_rate_limit_then_succeeds(monkeypatch):
    import urllib.error

    p = _provider(max_attempts=3)
    calls = {"n": 0}

    def flaky(payload):
        calls["n"] += 1
        if calls["n"] == 1:
            raise urllib.error.HTTPError("u", 429, "slow down", {}, None)
        return _reply(content="ok")

    monkeypatch.setattr(p, "_post", flaky)
    monkeypatch.setattr(P.time, "sleep", lambda *_: None)
    assert p.complete("s", "u") == "ok"
    assert P.STATS.transport_retries == 1
    assert P.STATS.failures == 0


def test_a_non_retryable_status_fails_immediately(monkeypatch):
    import urllib.error

    p = _provider(max_attempts=4)
    calls = {"n": 0}

    def not_found(payload):
        calls["n"] += 1
        raise urllib.error.HTTPError("u", 404, "no such model", {}, None)

    monkeypatch.setattr(p, "_post", not_found)
    monkeypatch.setattr(P.time, "sleep", lambda *_: None)
    with pytest.raises(RuntimeError, match="404"):
        p.complete("s", "u")
    assert calls["n"] == 1, "a 404 must not be retried"
    assert P.STATS.failures == 1


def test_exhausted_retries_raise_and_are_counted(monkeypatch):
    p = _provider(max_attempts=2)
    monkeypatch.setattr(p, "_post", lambda payload: (_ for _ in ()).throw(TimeoutError("t")))
    monkeypatch.setattr(P.time, "sleep", lambda *_: None)
    with pytest.raises(RuntimeError):
        p.complete("s", "u")
    assert P.STATS.failures == 1


def test_selection_defaults_by_model_family(monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "k")
    assert isinstance(P.get_provider("claude-opus-4-8"), P.AnthropicProvider)
    assert isinstance(P.get_provider("openai/gpt-oss-120b"), P.OpenAICompatProvider)


def test_explicit_provider_overrides_the_default(monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "k")
    p = P.get_provider("claude-opus-4-8", provider="ngc")
    assert isinstance(p, P.OpenAICompatProvider)
    assert p.base_url == P.NGC_BASE_URL


def test_ngc_without_a_key_fails_loudly():
    with pytest.raises(RuntimeError, match="NVIDIA_API_KEY"):
        P.get_provider("openai/gpt-oss-120b")


def test_unknown_provider_is_rejected():
    with pytest.raises(ValueError, match="unknown provider"):
        P.get_provider("some/model", provider="nope")


def test_stats_serialise():
    P.STATS.record("m", 5)
    d = P.STATS.as_dict()
    assert d["calls"] == 1 and d["total_tokens"] == 5 and d["by_model"] == {"m": 1}
    assert json.dumps(d)
