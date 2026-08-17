"""One code path for every proposer model in the Step 6 panel (PREREG_v4 section 7).

The panel compares model families, so the comparison is only about the models if
everything else is held identical: the same prompt, the same structured-output
validation, the same retry behaviour. That is why the providers sit behind one
interface and why nothing here is tunable per model. Section 7 forbids per-model
prompt tuning, and the temptation to do it is strongest on whichever model performs
worst, so there is deliberately nowhere to put it.

Two things the NGC-hosted models forced into the design.

First, they disagree about where the answer goes. `openai/gpt-oss-120b`,
`nvidia/nemotron-3-ultra-550b-a55b` and `thinkingmachines/inkling` populate a separate
`reasoning_content` field alongside `content`, and an early probe that read only
`content` recorded gpt-oss as returning an empty string when it had answered
perfectly well. A reader that misses this scores a working model as producing
nothing, and a parsing failure becomes indistinguishable from a model that cannot
follow the schema. Both fields are read.

Second, the free tier rate limits, so calls retry with backoff, and every call is
counted. `CallStats` is what section 7's proposal-validity and retry-rate columns are
computed from; a model that needed three attempts to emit valid JSON did not perform
as well as one that needed a single attempt, and the panel table has to show that.

Kept to the standard library on purpose. The OpenAI-compatible shape is a POST with a
JSON body, and adding an SDK would mean rebuilding the container image for no gain on
a host where disk is the binding constraint.
"""
from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field

NGC_BASE_URL = "https://integrate.api.nvidia.com/v1"
# A structure proposal is a 20,000-token structured generation, not a chat turn. The
# original 180 seconds was set against a trivial availability probe and timed out every
# call to the two largest panel models, which then read as "the model did not serve"
# when it was the client giving up. Sized for the slowest model that has to complete.
DEFAULT_TIMEOUT = 900
RETRY_STATUS = (408, 409, 425, 429, 500, 502, 503, 504)


@dataclass
class CallStats:
    """Per-model call accounting, for the panel's validity and retry columns."""

    calls: int = 0
    transport_retries: int = 0
    failures: int = 0
    total_tokens: int = 0
    empty_content_with_reasoning: int = 0
    by_model: dict[str, int] = field(default_factory=dict)

    def record(self, model: str, tokens: int = 0) -> None:
        self.calls += 1
        self.total_tokens += tokens
        self.by_model[model] = self.by_model.get(model, 0) + 1

    def as_dict(self) -> dict:
        return {"calls": self.calls, "transport_retries": self.transport_retries,
                "failures": self.failures, "total_tokens": self.total_tokens,
                "empty_content_with_reasoning": self.empty_content_with_reasoning,
                "by_model": dict(self.by_model)}


STATS = CallStats()


class Provider:
    """Anything that can turn a system and user prompt into text."""

    name = "provider"

    def complete(self, system: str, user: str, *, max_tokens: int = 20000,
                 effort: str = "high") -> str:
        raise NotImplementedError


class AnthropicProvider(Provider):
    """The existing path, unchanged, so the record stays comparable."""

    name = "anthropic"

    def __init__(self, model: str):
        self.model = model

    def complete(self, system: str, user: str, *, max_tokens: int = 20000,
                 effort: str = "high") -> str:
        from ..shared import llm
        text = llm.anthropic_reason(system, user, effort=effort, max_tokens=max_tokens)
        STATS.record(self.model)
        return text


class OpenAICompatProvider(Provider):
    """Any endpoint with the OpenAI chat-completions shape, including NGC."""

    name = "openai_compat"

    def __init__(self, model: str, base_url: str, api_key: str,
                 timeout: int = DEFAULT_TIMEOUT, max_attempts: int = 4):
        self.model = model
        self.base_url = base_url.rstrip("/")
        self.api_key = api_key
        self.timeout = timeout
        self.max_attempts = max_attempts

    def _post(self, payload: dict) -> dict:
        req = urllib.request.Request(
            f"{self.base_url}/chat/completions", data=json.dumps(payload).encode(),
            headers={"Authorization": f"Bearer {self.api_key}",
                     "Content-Type": "application/json"}, method="POST")
        with urllib.request.urlopen(req, timeout=self.timeout) as r:
            return json.loads(r.read().decode())

    def complete(self, system: str, user: str, *, max_tokens: int = 20000,
                 effort: str = "high") -> str:
        payload = {
            "model": self.model,
            "messages": [{"role": "system", "content": system},
                         {"role": "user", "content": user}],
            "max_tokens": max_tokens,
            "temperature": 0,
        }
        delay = 2.0
        last = ""
        for attempt in range(self.max_attempts):
            try:
                data = self._post(payload)
            except urllib.error.HTTPError as e:
                last = f"HTTP {e.code}: {e.read().decode('utf8', 'replace')[:160]}"
                if e.code not in RETRY_STATUS or attempt == self.max_attempts - 1:
                    STATS.failures += 1
                    raise RuntimeError(f"{self.model}: {last}") from e
            except (urllib.error.URLError, TimeoutError, OSError) as e:
                last = f"{type(e).__name__}: {e}"
                if attempt == self.max_attempts - 1:
                    STATS.failures += 1
                    raise RuntimeError(f"{self.model}: {last}") from e
            else:
                STATS.record(self.model, (data.get("usage") or {}).get("total_tokens", 0))
                return self._read(data)
            STATS.transport_retries += 1
            time.sleep(delay)
            delay *= 2
        STATS.failures += 1
        raise RuntimeError(f"{self.model}: {last}")

    def _read(self, data: dict) -> str:
        """Take the answer from wherever the model put it.

        Some models return the whole answer in `reasoning_content` and leave `content`
        empty. Both are returned, reasoning first, so the JSON block extractor
        downstream finds the structure regardless of which field carried it.
        """
        msg = ((data.get("choices") or [{}])[0].get("message")) or {}
        content = (msg.get("content") or "").strip()
        reasoning = (msg.get("reasoning_content") or "").strip()
        if reasoning and not content:
            STATS.empty_content_with_reasoning += 1
        if reasoning and content:
            return f"{reasoning}\n\n{content}"
        return content or reasoning


def _ngc_key() -> str:
    for var in ("NVIDIA_API_KEY", "NGC_API_KEY"):
        if os.environ.get(var):
            return os.environ[var]
    raise RuntimeError("Set NVIDIA_API_KEY for the NGC-hosted panel models.")


def get_provider(model: str | None = None, provider: str | None = None) -> Provider:
    """Build the provider for a model.

    Selection is explicit rather than guessed from the model string, because a wrong
    guess would silently route a panel model through the Anthropic client and report
    its failure as a model property.
    """
    model = model or os.environ.get("MMC_MODEL", "claude-opus-4-8")
    provider = (provider or os.environ.get("MMC_PROVIDER", "")).lower()
    if not provider:
        provider = "anthropic" if model.startswith("claude") else "ngc"
    if provider == "anthropic":
        return AnthropicProvider(model)
    if provider in ("ngc", "nvidia"):
        return OpenAICompatProvider(model, os.environ.get("MMC_BASE_URL", NGC_BASE_URL),
                                    _ngc_key())
    if provider in ("openai_compat", "openai"):
        base = os.environ.get("MMC_BASE_URL")
        key = os.environ.get("OPENAI_API_KEY") or _ngc_key()
        if not base:
            raise RuntimeError("Set MMC_BASE_URL for an OpenAI-compatible endpoint.")
        return OpenAICompatProvider(model, base, key)
    raise ValueError(f"unknown provider {provider!r}")
