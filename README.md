# Hypra

**Native Python/SymPy hypergeometric epsilon expansions, inspired by HypExp 2.**

No Mathematica installation is required. Hypra is an academic **partial native
port**, not yet a full replacement for HypExp/HPL. It computes exact symbolic
coefficients using recurrence relations and preserves unevaluated HPL words.
Unsupported cases fail explicitly; there is no silent external fallback.

## What works now

For rational slopes and real `0 < z < 1`, the native engine supports these
families to any requested finite order within the configured complexity limit:

- `2F1(1, b*eps; 1+c*eps; z)`: an all-order differential-equation/HPL recursion.
  This includes both the requested **P and Q** cases, not just their first terms.
- `2F1(1, 1/2+b*eps; 3/2+b*eps; z)`: an exact odd-polylogarithm expansion.
- `2F1(a*eps, 1/2+a*eps; 1/2; z)`: an exact even-binomial/logarithm expansion.

The complete expression may include Gamma functions, rational singular
prefactors, products and sums. Pole-aware order planning expands each
hypergeometric factor deeply enough before final Laurent truncation. Exact HPL
serialization, differentiation, selected elementary reductions, cache,
provenance, CLI and LaTeX output are included.

## Start

```sh
cd ~/project/hypra
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

- General `pFq`, general integer parameter shifts, most half-integer classes,
  full HPL reduction/analytic continuation, and generic integration remain to
  be ported. See [roadmap](docs/roadmap.md).
- Hypergeometric arguments must be regulator-independent and in the documented
  real domain. `z=0` is handled directly. Endpoint `z=1`, complex continuation,
  unknown valuations and essential/fractional-power regulator singularities
  fail closed.
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
