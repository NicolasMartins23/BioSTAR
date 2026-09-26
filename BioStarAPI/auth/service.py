from __future__ import annotations

import hashlib
import secrets
from datetime import UTC, datetime
from typing import Any

from fastapi import HTTPException, Request
from sqlalchemy import text

from BioStarAPI.database.connection import engine
from BioStarAPI.resources.exceptions import BioStarAPIError
from BioStarAPI.resources.messages import MessageCode

API_KEY_PREFIX = "bst_live_"
UNAUTH_DAILY_LIMIT = 8_640
AUTH_DAILY_LIMIT = 86_400
UNAUTH_IP_INTERVAL_SECONDS = 10
AUTH_IP_INTERVAL_SECONDS = 2


def _hash_api_key(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def generate_api_key(name: str | None = None) -> tuple[int, str]:
    raw_key = API_KEY_PREFIX + secrets.token_urlsafe(32)
    key_hash = _hash_api_key(raw_key)
    with engine.begin() as connection:
        row = connection.execute(
            text(
                """
                INSERT INTO api_keys (name, key_hash, key_prefix)
                VALUES (:name, :key_hash, :key_prefix)
                RETURNING id
                """
            ),
            {"name": name, "key_hash": key_hash, "key_prefix": raw_key[:16]},
        ).one()
    return int(row.id), raw_key


def revoke_api_key(key_id: int) -> bool:
    with engine.begin() as connection:
        result = connection.execute(
            text("UPDATE api_keys SET revoked_at = CURRENT_TIMESTAMP WHERE id = :id AND revoked_at IS NULL"),
            {"id": key_id},
        )
    return result.rowcount == 1


def resolve_api_key(raw_key: str | None) -> int | None:
    if not raw_key:
        return None
    with engine.connect() as connection:
        row = connection.execute(
            text(
                """
                SELECT id
                FROM api_keys
                WHERE key_hash = :key_hash
                  AND revoked_at IS NULL
                """
            ),
            {"key_hash": _hash_api_key(raw_key)},
        ).first()
    return int(row.id) if row else None


def _client_ip(request: Request) -> str:
    return request.client.host if request.client else "unknown"


def _check_short_rate_limit(request: Request, authenticated: bool) -> None:
    # One-process limiter is appropriate for the current single-container deployment.
    # Daily quotas are persisted in PostgreSQL and therefore survive restarts.
    import time

    interval = AUTH_IP_INTERVAL_SECONDS if authenticated else UNAUTH_IP_INTERVAL_SECONDS
    state = request.app.state.rate_limit_state
    key = ("auth" if authenticated else "anonymous", _client_ip(request))
    now = time.monotonic()
    previous = state.get(key)
    if previous is not None and now - previous < interval:
        retry_after = max(1, int(interval - (now - previous) + 0.999))
        raise BioStarAPIError(
            429,
            MessageCode.RATE_LIMIT_EXCEEDED,
            headers={"Retry-After": str(retry_after)},
        )
    state[key] = now


def _check_daily_limit(scope: str, scope_id: str | int, limit: int) -> None:
    today = datetime.now(UTC).date()
    with engine.begin() as connection:
        result = connection.execute(
            text(
                """
                INSERT INTO api_daily_usage (usage_date, scope, scope_id, request_count)
                VALUES (:usage_date, :scope, :scope_id, 1)
                ON CONFLICT (usage_date, scope, scope_id)
                DO UPDATE SET request_count = api_daily_usage.request_count + 1
                WHERE api_daily_usage.request_count < :limit
                RETURNING request_count
                """
            ),
            {
                "usage_date": today,
                "scope": scope,
                "scope_id": str(scope_id),
                "limit": limit,
            },
        ).first()
    if result is None:
        raise BioStarAPIError(429, MessageCode.DAILY_LIMIT_EXCEEDED)


def enforce_request_limits(request: Request, api_key_id: int | None) -> None:
    authenticated = api_key_id is not None
    _check_short_rate_limit(request, authenticated)
    if authenticated:
        _check_daily_limit("api_key", api_key_id, AUTH_DAILY_LIMIT)
    else:
        _check_daily_limit("server", "global", UNAUTH_DAILY_LIMIT)


def require_api_key(request: Request) -> int:
    api_key_id = resolve_api_key(request.headers.get("X-API-Key"))
    if api_key_id is None:
        raise BioStarAPIError(401, MessageCode.AUTHENTICATION_REQUIRED)
    return api_key_id
