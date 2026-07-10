"""Provenance receipt, gated value, and HMAC signing.

The receipt records the store computation behind every atlas-derived number that
enters the loop, so no fitted or observed value is untraceable. The gate enforces
the actor and tool boundary: the model sets structure, and every magnitude carries
a receipt from the fitter or the store.
"""
from __future__ import annotations

import hashlib
import hmac
import os
from datetime import UTC, datetime

from pydantic import BaseModel, Field


class Receipt(BaseModel):
    value: float | str | None
    source: str            # 'Zhu2025' | 'fit' | 'simulate' | 'Moonen2026'
    computation: str       # e.g. 'log_fc of GATA3 knockdown on IL13, Stim8hr'
    query: str
    retrieved_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    sig: str | None = None


class Gated(BaseModel):
    value: float | str | None
    receipt: Receipt


NO_RECEIPT = object()  # sentinel -> forces INSUFFICIENT / undetermined downstream


def sign(payload: bytes, key: bytes) -> str:
    return hmac.new(key, payload, hashlib.sha256).hexdigest()


def _key() -> bytes:
    k = os.environ.get("MMC_HMAC_KEY", "")
    return bytes.fromhex(k) if k else b"dev-key-not-for-production"
