# Hypra

**Native Python/SymPy hypergeometric epsilon expansions, inspired by HypExp 2.**

Native HPL functionality is also informed by HPL 2.0.
No Mathematica installation is required. Hypra is an academic **partial native
implementation**, not yet a full replacement for HypExp/HPL. It computes exact symbolic
coefficients using recurrence relations and preserves unevaluated HPL words.
Unsupported cases fail explicitly; there is no silent external fallback.

## What works now

The implementation now has general-purpose native components, rather than a
list of hard-coded low-order answers:

- Integer-base `pF(p-1)` with exact rational or guarded symbolic polynomial slopes: strict nested sums,
  contiguous parameter shifts, and finite-prefix/tail reindexing for moving
  negative upper bases. Moving singular lower parameters retain their poles.
- Matched-half `pF(p-1)`: equal counts of half-integer upper/lower parameters,
  colored strict sums and HPLs at `sqrt(z)`, with integer-offset reductions.
- Four canonical half-integer differential patterns, including unmatched
  and all-half families. Combined with the sum kernels, their canonical
  parameter-count coverage matches the HypExp 2 whitelist; this is not a
  claim of complete upstream API or analytic-continuation parity.
- Exact fixed terminating sums and convergent Gauss `2F1` unit-endpoint
  rewriting, with explicit rejection when convergence cannot be established.
- HPL shuffle products, rational-HPL integration with a verified finite origin
  constant, exact local HPL series, shuffle logarithm extraction, convergent
  binary endpoints/MZVs, and branch-explicit negation/reflection transforms.

The main nonterminating argument domain is real `0 < z < 1`; use exact
affine exact parameters; singular symbolic lower slopes need explicit nonzero assumptions. See the [capability map](docs/scope.md) for module contracts,
public-dispatch distinctions and restrictions. The P/Q recurrence and two
special half-integer closed forms remain useful optimized cases, not the
extent of the native algorithms.

The complete expression may include Gamma functions, rational singular
prefactors, products and sums. Pole-aware order planning expands each
hypergeometric factor deeply enough before final Laurent truncation. Exact HPL
serialization, differentiation, selected elementary reductions, cache,
provenance, CLI and LaTeX output are included.

## Start

```sh
cd ~/projects/hypra
uv sync
uv run hypra backend-info
uv run python scripts/validate.py --output-dir reports/local
```

Python 3.12+ is required. Only Python dependencies are installed.

```python
import sympy as s
from hypra import expand_epsilon

eps, w = s.symbols("eps w")
P = s.hyper([1, -eps], [1 - 2*eps], w)
result = expand_epsilon(P, eps, through=2, assumptions=s.And(w > 0, w < 1))
print(result.expression)
# 1 + eps*log(1-w)
#   + eps**2*(-2*polylog(2,w)-log(1-w)**2/2) + O(eps**3)
print(result.latex())
```

`through=k` is inclusive and includes every Laurent pole term. For example,
`expand_epsilon(gamma(eps)*P, eps, through=2, assumptions=...)` requests P one
order deeper, through cubic order. `result.polynomial` excludes the remainder;
`result.coefficient(-1)` gives the pole coefficient; `result.to_dict()` preserves
exact ASTs and reproducibility metadata. `sympy.sympify(result)` includes the
remainder.

## HPL support

```python
from hypra import HPL, encode_expr, decode_expr
from hypra.symbolic import reduce_known_hpl
h = HPL((0, 1, -1), w)
assert decode_expr(encode_expr(h)) == h
assert HPL.from_compressed((2, -1), w) == h
print(s.diff(h, w))
print(s.latex(h))
```

Expanded letters `{-1,0,1}` use kernels `1/(1+t)`, `1/t`, `1/(1-t)`;
`H_1(z)=-log(1-z)`. Construction/serialization preserves words. Known-word
reduction is explicit and does not claim arbitrary analytic continuation.

Additional native HPL operations include exact at-i values for low weights,
explicit upper/lower real-cut numerical values, ordinary-point Taylor series,
classical reduction of all cubic words with at most two distinct letters,
alternating endpoints, MZV/EulerSum stuffle products and explicitly named
shuffle-regularization APIs. Each has documented domain/resource limits;
regularized constants are never silently treated as ordinary finite limits.
See the [HPL 2 capability map](docs/parity-map.md) for remaining gaps.

## Validate and produce LaTeX

```sh
uv run pytest -q
uv run ruff check .
uv run python scripts/validate.py --numeric --output-dir reports/local
uv run hypra expand examples/gamma_q_request.json \
  --output result.json --latex result.tex
```

The native validation **runs the Python engine** on P, Q and the complete J4
integrand, then uses five independently derived logarithmic moments to integrate
its coefficients. The mathematical documentation proves the endpoint exchange
for J4 specifically. This yields

\[
J_4=2\zeta(3)\epsilon^2+\frac{11\pi^4}{60}\epsilon^3+O(\epsilon^4).
\]

This is not a general symbolic integration engine. No expected coefficient is
injected into the native algorithm. Numerical quadrature is secondary evidence.

See [mathematics](docs/mathematics.md), its [LaTeX source](docs/mathematics.tex)
and compiled [PDF](docs/mathematics.pdf). Runnable examples include
`examples/gamma_q_request.json` and its actual native JSON/LaTeX outputs.

## Boundaries

- Nonterminating general `pFq` outside `pF(p-1)`, uncovered half-integer
  patterns, full HPL reduction/analytic continuation, and generic integration
  remain unsupported. See [roadmap](docs/roadmap.md).
- Hypergeometric arguments must be regulator-independent and in the documented
  real domain. `z=0` is handled directly. Only documented terminating/Gauss
  identities cover `z=1`; arbitrary unit endpoints and complex continuation
  are not supported. Unknown valuations and essential/fractional-power
  regulator singularities fail closed.
- Use exact `Rational` values, never binary floats. Unknown functions, invalid
  identifiers and contradictory/unsupported domain assumptions are rejected.
- `NativeConfig(max_internal_order=12)` is the default complexity limit; it can
  be raised explicitly. HPL word counts can grow exponentially. There is no
  claim that arbitrary large requests will be fast.
- The optional external Wolfram comparison bridge is **experimental, unverified
  and not needed**. Its strict remainder checks reject upstream paths returning
  bare polynomials; it is not advertised as a working comparison oracle.

## Upstream source and academic use

Official HypExp/HPL archives were retrieved into Git-ignored `reference/upstream`
and inspected read-only. [Source manifest and algorithm map](docs/upstream-source.md)
record exact origins/hashes and translation priorities. The current native
formulas are independently implemented mathematical algorithms with citations,
not copied upstream tables. The CPC licence restricts redistribution of original
or modified source; no upstream code is bundled or published. Review permission
before publishing a direct translation. No project license or registry/trademark
availability for the working name **Hypra** is claimed yet.

The NNLO project was read only and remains unchanged.

References: [HypExp](https://arxiv.org/abs/hep-ph/0507094),
[HypExp 2](https://arxiv.org/abs/0708.2443).

## HPL tools and bounded jobs

```python
from hypra import (HPL, hpl_at_one, extract_hpl_logs, negate_hpl, reflect_hpl,
                   NativeConfig)

extract_hpl_logs(HPL((1, 0), w))
hpl_at_one((0, 0, 1))  # zeta(3)
reflect_hpl((0, 1), w, domain="0<z<1")
negate_hpl((1, 0), w, domain="0<z<1")
config = NativeConfig(timeout_seconds=30)  # killable backend process
```

Transformations use explicit domains and reject unsupported branches. General
MZVs remain exact rather than being reduced with unverified tables. See
[HPL endpoints](docs/hpl-endpoints.md), [log extraction](docs/hpl-reduction.md),
and [argument transformations](docs/hpl-transforms.md).

For numerical HPL values, `hpl_evalf(word, z, dps=50)` defaults to the real
interval. Select `domain="principal-disk"` explicitly for complex points inside
the slit unit disk. `continue_hpl(word, z, domain="principal-cut-plane",
path="radial", dps=50)` continues numerically along the specified off-cut
radial path; it reports propagated truncation and precision-stability diagnostics.
These numerical diagnostics do not certify floating-point roundoff.

## Browse the documentation locally

```sh
uv sync
uv run mkdocs serve --dev-addr 127.0.0.1:8000
```

Open http://127.0.0.1:8000/ on the same computer. The documentation includes
locally bundled math rendering and works without external browser resources.
