"""Strict JSON loading shared by the CLI and both backends."""

import json
from typing import Any


def strict_json_loads(text: str) -> Any:
    """Reject duplicate object keys and JSON's nonstandard nonfinite constants."""

    def unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> Any:
        raise ValueError(f"Nonfinite JSON constant: {value}")

    return json.loads(text, object_pairs_hook=unique_pairs, parse_constant=reject_constant)
