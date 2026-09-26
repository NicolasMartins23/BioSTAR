from __future__ import annotations

from collections.abc import Mapping

from fastapi import Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from BioStarAPI.resources.messages import MessageCode, get_message
from BioStarAPI.resources.responses import APIResponse


class BioStarAPIError(Exception):
    def __init__(
        self,
        status_code: int,
        message_code: MessageCode,
        headers: dict[str, str] | None = None,
        **parameters: object,
    ) -> None:
        self.status_code = status_code
        self.headers = headers or {}
        self.message = get_message(message_code, **parameters)
        super().__init__(self.message.message)


def register_exception_handlers(app) -> None:
    @app.exception_handler(BioStarAPIError)
    async def handle_api_error(request: Request, exception: BioStarAPIError) -> JSONResponse:
        return JSONResponse(
            status_code=exception.status_code,
            headers=exception.headers,
            content=APIResponse[None](
                data=None,
                message=exception.message,
            ).model_dump(mode="json"),
        )

    @app.exception_handler(RequestValidationError)
    async def handle_validation_error(
        request: Request,
        exception: RequestValidationError,
    ) -> JSONResponse:
        return JSONResponse(
            status_code=422,
            content=APIResponse[None](
                data=None,
                message=get_message(MessageCode.VALIDATION_ERROR),
            ).model_dump(mode="json"),
        )
