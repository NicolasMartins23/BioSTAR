from __future__ import annotations

import os

from fastapi import APIRouter, Header
from pydantic import BaseModel, Field

from BioStarAPI.auth.service import generate_api_key, revoke_api_key
from BioStarAPI.resources.exceptions import BioStarAPIError
from BioStarAPI.resources.messages import MessageCode, get_message
from BioStarAPI.resources.responses import APIResponse

router = APIRouter(prefix="/api/auth", tags=["auth"])


class CreateAPIKeyRequest(BaseModel):
    name: str | None = Field(default=None, max_length=100)


class CreateAPIKeyResponse(BaseModel):
    id: int
    name: str | None
    api_key: str
    


@router.post("/keys", response_model=APIResponse[CreateAPIKeyResponse])
def create_api_key(
    request: CreateAPIKeyRequest,
    x_admin_key: str | None = Header(default=None, alias="X-Admin-Key"),
) -> CreateAPIKeyResponse:
    expected = os.getenv("BIOSTAR_AUTH_ADMIN_KEY")
    if not expected or not x_admin_key or not secrets_compare(x_admin_key, expected):
        raise BioStarAPIError(401, MessageCode.INVALID_ADMIN_KEY)
    key_id, api_key = generate_api_key(request.name)
    return APIResponse(
        data=CreateAPIKeyResponse(id=key_id, name=request.name, api_key=api_key),
        message=get_message(MessageCode.API_KEY_CREATED),
    )


@router.delete("/keys/{key_id}")
def delete_api_key(
    key_id: int,
    x_admin_key: str | None = Header(default=None, alias="X-Admin-Key"),
) -> dict[str, object]:
    expected = os.getenv("BIOSTAR_AUTH_ADMIN_KEY")
    if not expected or not x_admin_key or not secrets_compare(x_admin_key, expected):
        raise HTTPException(status_code=401, detail="Invalid authentication administrator key")
    if not revoke_api_key(key_id):
        raise BioStarAPIError(404, MessageCode.API_KEY_NOT_FOUND)
    return APIResponse(
        data={"revoked": True, "id": key_id},
        message=get_message(MessageCode.API_KEY_REVOKED),
    )


def secrets_compare(left: str, right: str) -> bool:
    import secrets

    return secrets.compare_digest(left, right)
