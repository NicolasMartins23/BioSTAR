from __future__ import annotations

from typing import Generic, TypeVar

from pydantic import BaseModel

from BioStarAPI.resources.messages import MessageResource

T = TypeVar("T")


class APIResponse(BaseModel, Generic[T]):
    data: T | None = None
    message: MessageResource | None = None
