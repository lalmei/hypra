"""Hypra's real Typer command-line application: exact JSON input, never eval."""

import json
import sys
from enum import Enum
from pathlib import Path
from typing import Annotated

import typer

from .api import backend_info, expand_epsilon
from .backend import BackendError, _strict_json_loads
from .errors import HypraError
from .serialization import SerializationError, decode_expr

app = typer.Typer(
    name="hypra",
    help="Native exact hypergeometric epsilon expansions from structured JSON AST input.",
    add_completion=False,
    no_args_is_help=False,
    pretty_exceptions_enable=False,
)


class Backend(str, Enum):
    """Select the native engine or the optional external comparison bridge."""

    native = "native"
    wolfram = "wolfram"


_ERRORS = (BackendError, HypraError, SerializationError, OSError, ValueError, TypeError, KeyError)


def _error(exc):
    print(json.dumps({"status": "error", "type": type(exc).__name__, "message": str(exc)}),
          file=sys.stderr)
    raise typer.Exit(code=2) from exc


def _service(ctx, name, default):
    # main(argv) can inject its imported functions, retaining backwards-compatible
    # programmatic testing without a second/fake CLI implementation.
    return (ctx.obj or {}).get(name, default)


@app.command("backend-info")
def backend_info_command(
    ctx: typer.Context,
    backend: Annotated[Backend, typer.Option(help="Backend whose availability and capabilities to report.")]
    = Backend.native,
):
    """Report backend capabilities as JSON; exit 2 when unavailable."""
    try:
        info = _service(ctx, "backend_info", backend_info)(backend=backend.value)
        print(json.dumps(info, indent=2))
        if not info["available"]:
            raise typer.Exit(code=2)
    except _ERRORS as exc:
        _error(exc)


@app.command("expand")
def expand_command(
    ctx: typer.Context,
    request: Annotated[Path, typer.Argument(help="JSON file with exact expression and epsilon ASTs.")],
    backend: Annotated[Backend, typer.Option(help="Expansion engine; no silent external fallback.")]
    = Backend.native,
    output: Annotated[Path | None, typer.Option(help="Write result JSON here instead of stdout.")]
    = None,
    latex: Annotated[Path | None, typer.Option(help="Also write a LaTeX display equation.")]
    = None,
):
    """Expand an exact JSON AST request through its inclusive Laurent order.

    Request keys: expression, epsilon, through (default 3), assumptions.
    Input is data, never Python or Mathematica source. Errors are JSON on
    stderr with exit status 2. Successful results retain exact provenance.
    """
    try:
        payload = _strict_json_loads(request.read_text())
        if not isinstance(payload, dict) or set(payload) - {
            "expression", "epsilon", "through", "assumptions",
        }:
            raise ValueError("Request must contain only expression, epsilon, through, assumptions")
        result = _service(ctx, "expand_epsilon", expand_epsilon)(
            decode_expr(payload["expression"]),
            decode_expr(payload["epsilon"]),
            through=payload.get("through", 3),
            backend=backend.value,
            assumptions=decode_expr(payload["assumptions"]) if "assumptions" in payload else None,
        )
        serialized = json.dumps(result.to_dict(), indent=2) + "\n"
        if output:
            output.write_text(serialized)
        else:
            print(serialized, end="")
        if latex:
            latex.write_text(
                "% Exact Hypra result; see JSON for provenance.\n\\[\n" + result.latex() + "\n\\]\n"
            )
    except _ERRORS as exc:
        _error(exc)
