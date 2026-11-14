# Command-line reference

Generated from the actual Typer application used by `hypra` and `python -m hypra`.
Input expressions use exact JSON AST, never evaluated Python source strings.
The CLI uses the default native argument and unit-endpoint policies. Use
[`expand_epsilon` with `NativeConfig`](python-api.md) to select a principal
complex domain or meromorphic parameter continuation.

::: mkdocs-typer2
    :module: hypra.cli
    :name: hypra
