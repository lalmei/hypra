"""Assemble the Typer application and handle process exit semantics."""

from collections.abc import Sequence

import typer

try:
    # Recent Typer releases vendor Click; older supported releases use Click.
    from typer._click.exceptions import Abort, ClickException
except ImportError:
    from click import Abort, ClickException

from ..api import backend_info, expand_epsilon
from .commands import Commands


class CliApplication:
    """Register Hypra's commands with injectable expansion and status services."""

    def __init__(self, *, get_backend_info=backend_info, expand=expand_epsilon):
        self.app = typer.Typer(
            name="hypra",
            help="Native exact hypergeometric epsilon expansions from structured JSON AST input.",
            add_completion=False,
            no_args_is_help=False,
            pretty_exceptions_enable=False,
        )
        commands = Commands(get_backend_info, expand)
        self.app.command("backend-info")(commands.backend_info)
        self.app.command("expand")(commands.expand)

    def run(self, argv: Sequence[str] | None = None) -> int:
        """Return the command status; usage errors and aborts raise SystemExit."""
        try:
            result = self.app(args=argv, prog_name="hypra", standalone_mode=False)
            return result if isinstance(result, int) else 0
        except ClickException as exc:
            exc.show()
            raise SystemExit(exc.exit_code) from exc
        except Abort as exc:
            raise SystemExit(1) from exc


_cli = CliApplication()
app = _cli.app


def main(argv: Sequence[str] | None = None) -> int:
    """Invoke the same application used by the console script and documentation."""
    return _cli.run(argv)
