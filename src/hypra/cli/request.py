"""Read data-only expansion requests from disk."""

from dataclasses import dataclass
from pathlib import Path
from typing import Self

from sympy import Basic

from .._json import strict_json_loads
from ..serialization import decode_expr


@dataclass(frozen=True)
class ExpansionRequest:
    expression: Basic
    epsilon: Basic
    through: int
    assumptions: Basic | None

    @classmethod
    def read(cls, path: Path) -> Self:
        payload = strict_json_loads(path.read_text())
        if not isinstance(payload, dict) or set(payload) - {
            "expression", "epsilon", "through", "assumptions",
        }:
            raise ValueError("Request must contain only expression, epsilon, through, assumptions")
        return cls(
            expression=decode_expr(payload["expression"]),
            epsilon=decode_expr(payload["epsilon"]),
            through=payload.get("through", 3),
            assumptions=decode_expr(payload["assumptions"]) if "assumptions" in payload else None,
        )
