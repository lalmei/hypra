# Hypra

## Exact epsilon expansions, natively in Python

Hypra expands hypergeometric expressions with SymPy, retaining exact rational,
Gamma, zeta and harmonic-polylogarithm coefficients. Mathematica is **not**
required for the native engine.

!!! warning "Research software: partial native implementation"
    Hypra is inspired by HypExp/HypExp 2, not yet a complete replacement.
    Supported cases have explicit domain and resource checks. Unsupported
    cases fail rather than silently returning incomplete expansions.

### Start exploring

- **[Run an expansion](api-examples.md)**: inclusive Laurent orders, singular
  prefactors, exact result objects and LaTeX output.
- **[Understand the coverage](scope.md)**: native kernels, parameter families,
  public-dispatch restrictions and branch conventions.
- **[Read the derivation](mathematics.md)** or open the
  **[typeset mathematical report](mathematics.pdf)**.
- **[Evaluate HPLs](hpl-numerics.md)** with explicit real or principal-disk
  domains, and inspect numerical accuracy diagnostics.

### Current implementation

| Area | Status |
| --- | --- |
| Integer-base `pF(p-1)` and supported half-integer families | Native exact algorithms, with documented limits |
| Laurent poles and Gamma prefactors | Pole-aware internal order planning |
| HPL algebra, endpoints and argument transformations | Exact operations with explicit scope and branch rules |
| HPL numerical evaluation and radial continuation | Approximate high precision; no interval roundoff certification |
| Optional Wolfram comparison bridge | Experimental and not validated as a real-backend oracle |
| Full HypExp/HPL parity | Not claimed; see the [roadmap](roadmap.md) |

The primary nonterminating hypergeometric domain is real `0 < z < 1`.
Separate HPL tools can expose broader domains only through explicit options.
Use exact rational parameter data; numerical validation is secondary evidence,
not a substitute for symbolic derivations.

### Local documentation

This site runs only on this Mac at `http://127.0.0.1:8000/`. Search, styles,
fonts and mathematical rendering are local; no account or hosted site is used.

To restart it from the project directory:

```sh
uv run mkdocs serve --dev-addr 127.0.0.1:8000
```

To build static HTML: `uv run mkdocs build`. Generated files are not package
source. Edits to these Markdown pages refresh the development server.
