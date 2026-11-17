# Hypra

## Exact epsilon expansions, natively in Python

Hypra expands hypergeometric expressions with SymPy, retaining exact rational,
Gamma, zeta and harmonic-polylogarithm coefficients. Mathematica is **not**
required for the native engine.

!!! note "Research software: explicit mathematical contracts"
    Hypra implements the operations audited against HypExp 2 and HPL 2.
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
| [Unit endpoints](endpoints.md) | Convergent Gauss and supported higher-rank values, with opt-in meromorphic parameter continuation |
| [Definite logarithmic integrals](python-api.md#definite-logarithmic-integrals) | Exact beta-log moments and real or explicit principal nonreal mixed-log moments |
| HPL algebra, endpoints and argument transformations | Exact operations with explicit scope and branch rules |
| HPL numerical evaluation and radial continuation | Approximate high precision; no interval roundoff certification |
| Optional Wolfram comparison bridge | Live regression coverage against pinned HypExp 2 and HPL 2 |
| HypExp/HPL operation coverage | See the [port audit](port-audit.md), including symbolic normal-form distinctions |

The default nonterminating interior domain is real `0 < z < 1`.
Select a [principal complex domain](scope.md#principal-complex-arguments)
explicitly for exact numerical hypergeometric arguments. The HPL numerical
APIs have their own domain and cut-side options.

### Local documentation

Serve the documentation at `http://127.0.0.1:8000/` from the project
directory. Search, styles, fonts and math rendering use bundled assets.

```sh
uv run properdocs serve --dev-addr 127.0.0.1:8000
```

To build static HTML: `uv run properdocs build`. Generated files are not package
source. Edits to these Markdown pages refresh the development server.
