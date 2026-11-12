"""Run the real Typer CLI, preserving the programmatic main(argv) entrypoint."""

import click

from .api import backend_info, expand_epsilon
from .cli import app


def main(argv=None):
    """Invoke Hypra's CLI and return its execution exit status.

    As with the previous argparse entrypoint, malformed command syntax raises
    SystemExit; handled request/backend errors return 2 and print JSON stderr.
    """
    try:
        result = app(
            args=argv,
            prog_name="hypra",
            standalone_mode=False,
            obj={"backend_info": backend_info, "expand_epsilon": expand_epsilon},
        )
        return result if isinstance(result, int) else 0
    except click.ClickException as exc:
        exc.show()
        raise SystemExit(exc.exit_code) from exc
    except click.Abort as exc:
        raise SystemExit(1) from exc


if __name__ == "__main__":
    raise SystemExit(main())
