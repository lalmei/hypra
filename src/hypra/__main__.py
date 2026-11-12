"""Command line interface using structured exact input, never eval."""

import argparse
import json
import sys
from pathlib import Path

from .api import backend_info, expand_epsilon
from .backend import BackendError, _strict_json_loads
from .errors import HypraError
from .serialization import SerializationError, decode_expr


def main(argv=None):
    parser = argparse.ArgumentParser(prog="hypra", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    info_parser = commands.add_parser("backend-info", help="Report backend capabilities")
    info_parser.add_argument("--backend", choices=("native", "wolfram"), default="native")
    expand = commands.add_parser("expand", help="Expand an exact JSON AST request")
    expand.add_argument("request", type=Path)
    expand.add_argument("--backend", choices=("native", "wolfram"), default="native")
    expand.add_argument("--output", type=Path, help="Write result JSON instead of stdout")
    expand.add_argument("--latex", type=Path, help="Also write a LaTeX display equation")
    args = parser.parse_args(argv)
    try:
        if args.command == "backend-info":
            info = backend_info(backend=args.backend)
            print(json.dumps(info, indent=2))
            return 0 if info["available"] else 2
        request = _strict_json_loads(args.request.read_text())
        if not isinstance(request, dict) or set(request) - {
            "expression",
            "epsilon",
            "through",
            "assumptions",
        }:
            raise ValueError("Request must contain only expression, epsilon, through, assumptions")
        result = expand_epsilon(
            decode_expr(request["expression"]),
            decode_expr(request["epsilon"]),
            through=request.get("through", 3),
            backend=args.backend,
            assumptions=decode_expr(request["assumptions"]) if "assumptions" in request else None,
        )
        output = json.dumps(result.to_dict(), indent=2) + "\n"
        if args.output:
            args.output.write_text(output)
        else:
            print(output, end="")
        if args.latex:
            args.latex.write_text(
                "% Exact Hypra result; see JSON for provenance.\n\\[\n" + result.latex() + "\n\\]\n"
            )
        return 0
    except (
        BackendError,
        HypraError,
        SerializationError,
        OSError,
        ValueError,
        TypeError,
        KeyError,
    ) as exc:
        print(
            json.dumps({"status": "error", "type": type(exc).__name__, "message": str(exc)}),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
