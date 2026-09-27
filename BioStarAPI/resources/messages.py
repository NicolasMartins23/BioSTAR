from __future__ import annotations

from enum import Enum

from pydantic import BaseModel


class MessageCode(str, Enum):
    API_KEY_CREATED = "api_key_created"
    API_KEY_REVOKED = "api_key_revoked"
    API_KEY_NOT_FOUND = "api_key_not_found"
    AUTHENTICATION_REQUIRED = "authentication_required"
    INVALID_API_KEY = "invalid_api_key"
    INVALID_ADMIN_KEY = "invalid_admin_key"

    INVALID_FASTA = "invalid_fasta"
    FASTA_SINGLE_SEQUENCE_REQUIRED = "fasta_single_sequence_required"
    EMPTY_SEQUENCE = "empty_sequence"
    SEQUENCE_TOO_LONG = "sequence_too_long"
    INVALID_SEQUENCE = "invalid_sequence"
    MUTATION_SEQUENCES_SAME_LENGTH = "mutation_sequences_same_length"
    MUTATION_SEQUENCES_MULTIPLE_OF_THREE = "mutation_sequences_multiple_of_three"
    VALIDATION_ERROR = "validation_error"

    RATE_LIMIT_EXCEEDED = "rate_limit_exceeded"
    DAILY_LIMIT_EXCEEDED = "daily_limit_exceeded"

    STANDARD_GENETIC_CODE_NOT_SEEDED = "standard_genetic_code_not_seeded"
    AMINO_ACID_DATA_NOT_SEEDED = "amino_acid_data_not_seeded"


class MessageResource(BaseModel):
    code: MessageCode
    message: str
    details: list[dict[str, object]] | None = None


_MESSAGES: dict[MessageCode, str] = {
    MessageCode.API_KEY_CREATED: "API key created successfully.",
    MessageCode.API_KEY_REVOKED: "API key revoked successfully.",
    MessageCode.API_KEY_NOT_FOUND: "API key not found.",
    MessageCode.AUTHENTICATION_REQUIRED: "A valid X-API-Key is required.",
    MessageCode.INVALID_API_KEY: "The provided API key is invalid or revoked.",
    MessageCode.INVALID_ADMIN_KEY: "Invalid authentication administrator key.",
    MessageCode.INVALID_FASTA: "Invalid FASTA sequence.",
    MessageCode.FASTA_SINGLE_SEQUENCE_REQUIRED: "Exactly one FASTA sequence is required.",
    MessageCode.EMPTY_SEQUENCE: "{kind} sequence cannot be empty.",
    MessageCode.SEQUENCE_TOO_LONG: "{kind} sequence cannot exceed {max_length} characters.",
    MessageCode.INVALID_SEQUENCE: "Invalid {kind} sequence characters: {characters}.",
    MessageCode.MUTATION_SEQUENCES_SAME_LENGTH: "Reference and sequence must have the same length.",
    MessageCode.MUTATION_SEQUENCES_MULTIPLE_OF_THREE: "Reference and sequence lengths must be multiples of 3.",
    MessageCode.VALIDATION_ERROR: "The request contains invalid or missing fields.",
    MessageCode.RATE_LIMIT_EXCEEDED: "Rate limit exceeded.",
    MessageCode.DAILY_LIMIT_EXCEEDED: "Daily request limit exceeded.",
    MessageCode.STANDARD_GENETIC_CODE_NOT_SEEDED: "Standard genetic code is not seeded.",
    MessageCode.AMINO_ACID_DATA_NOT_SEEDED: "Amino acid reference data is not seeded.",
}


def get_message(code: MessageCode, **parameters: object) -> MessageResource:
    template = _MESSAGES[code]
    return MessageResource(code=code, message=template.format(**parameters))
