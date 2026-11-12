# Building and maintaining the documentation

The theme is **MaterialX**, installed as `mkdocs-materialx` from PyPI and
configured as `theme.name: materialx`. It is the project at
https://github.com/jaywhj/mkdocs-materialx, not the similarly named original
`mkdocs-material` distribution. Do not install both into the same environment:
they provide overlapping Python/theme files.

```sh
uv sync --locked --group docs
uv run --no-sync mkdocs build --strict
uv run --no-sync python scripts/check_docs_build.py
uv run --no-sync mkdocs serve --dev-addr 127.0.0.1:8000
```

The `docs` dependency group is optional for mathematical package use. The lock
file fixes the versions used by CI. The repository includes both `dev` and `docs` in its default uv groups, so
ordinary development commands do not remove the running site's plugins.
Use `--no-group docs` only when intentionally working without documentation tools. Do not start a second development server on the same port if the
existing managed local server is already running.

## Documentation components

- **KaTeX:** locally vendored CSS, JavaScript and fonts. `arithmatex` preserves
  LaTeX mathematics for browser rendering. Both inline `\(...\)` and display
  `\[...\]` forms are supported; no remote math-rendering service is used.
- **mkdocstrings-python:** collects actual Python signatures and docstrings for
  the [API reference](python-api.md), rather than maintaining a second signature
  list manually.
- **mkdocs-typer2:** documents the real `hypra.cli:app` Typer application. The
  CLI was migrated from argparse; there is no fake documentation-only Typer
  wrapper. The native renderer and Typer must be upgraded together and tested:
  the initial combination with Typer 0.29.1 failed in option introspection.
- **mermaid2:** renders architecture diagrams using a pinned local Mermaid ESM
  bundle. Its dependency graph, source hash and MIT license are under
  `assets/mermaid/`. Standard diagram rendering does not depend on a CDN;
  optional externally supplied ELK worker functionality is not bundled.
- **LaTeX:** the mathematical report remains available as both
  [editable TeX](mathematics.tex) and [PDF](mathematics.pdf). KaTeX handles web
  equations; it does not compile the PDF report.

## CI is build-only

`.github/workflows/docs.yml` installs locked dependencies, performs a strict
build, checks generated API/CLI pages and local assets, and uploads the static
site as a short-lived Actions artifact. It has read-only repository permissions.
It does **not** enable GitHub Pages, deploy publicly, push a `gh-pages` branch,
or change repository settings.

`check_docs_build.py` checks generated structure and local dependency files.
It is not a browser-layout or JavaScript execution test. Browser smoke checks
should additionally confirm rendered `.katex` elements and Mermaid SVGs when
changing the theme or frontend libraries. Keep a mathematical fallback visible
if rendering fails: the TeX source and downloadable PDF remain available.
