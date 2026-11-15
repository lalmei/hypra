"""Typer commands for backend discovery and exact expansion."""

import json
import sys
from enum import Enum
from pathlib import Path
from typing import Annotated, NoReturn

import typer

from ..backend import BackendError
from ..errors import HypraError
from .request import ExpansionRequest


class Backend(str, Enum):
    """Select the native engine or the optional external comparison bridge."""

    native = "native"
    wolfram = "wolfram"


_ERRORS = (BackendError, HypraError, OSError, ValueError, TypeError, KeyError)


def _error(exc: Exception) -> NoReturn:
    print(
        json.dumps({"status": "error", "type": type(exc).__name__, "message": str(exc)}),
        file=sys.stderr,
    )
    raise typer.Exit(code=2) from exc


class Commands:
    """Translate command inputs and outputs around Hypra's Python API."""

    def __init__(self, get_backend_info, expand):
        self._get_backend_info = get_backend_info
        self._expand = expand

    def backend_info(
        self,
        backend: Annotated[
            Backend, typer.Option(help="Backend whose availability and capabilities to report.")
        ] = Backend.native,
    ) -> None:
        """Report backend capabilities as JSON; exit 2 when unavailable."""
        try:
            info = self._get_backend_info(backend=backend.value)
            print(json.dumps(info, indent=2))
            if not info["available"]:
                raise typer.Exit(code=2)
        except _ERRORS as exc:
            _error(exc)

    def expand(
        self,
        request: Annotated[
            Path, typer.Argument(help="JSON file with exact expression and epsilon ASTs.")
        ],
        backend: Annotated[
            Backend, typer.Option(help="Expansion engine; no silent external fallback.")
        ] = Backend.native,
        output: Annotated[
            Path | None, typer.Option(help="Write result JSON here instead of stdout.")
        ] = None,
        latex: Annotated[
            Path | None, typer.Option(help="Also write a LaTeX display equation.")
        ] = None,
    ) -> None:
        """Expand an exact JSON AST request through its inclusive Laurent order.

        Request keys: expression, epsilon, through (default 3), assumptions.
        Input is data, never Python or Mathematica source. Errors are JSON on
        stderr with exit status 2. Successful results retain exact provenance.
        """
        try:
            payload = ExpansionRequest.read(request)
            result = self._expand(
                payload.expression,
                payload.epsilon,
                through=payload.through,
                backend=backend.value,
                assumptions=payload.assumptions,
            )
            serialized = json.dumps(result.to_dict(), indent=2) + "\n"
            if output:
                output.write_text(serialized)
            else:
                print(serialized, end="")
            if latex:
                latex.write_text(
                    f"% Exact Hypra result; see JSON for provenance.\n\\[\n{result.latex()}\n\\]\n"
                )
        except _ERRORS as exc:
            _error(exc)
