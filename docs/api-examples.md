# API examples

## Exact epsilon expansion

```python
import sympy as s
from hypra import expand_epsilon

eps, z = s.symbols("eps z")
P = s.hyper([1, -eps], [1 - 2*eps], z)
result = expand_epsilon(
    P, eps, through=2, assumptions=s.And(z > 0, z < 1)
)
print(result.expression)       # Includes O(eps**3)
print(result.polynomial)       # No remainder
print(result.coefficient(2))   # Exact quadratic coefficient
print(result.latex())
```

`through=2` includes **every Laurent term through epsilon squared**, including
negative powers. For `gamma(eps)*P`, the planner requests a deeper expansion of
`P` before multiplying and truncating; do not expand factors independently at
insufficient depth.

## Exact HPL objects and endpoint constants

```python
from hypra import HPL, encode_expr, decode_expr
from hypra import hpl_at_one, hpl_alternating_at_one

h = HPL((0, 1, -1), z)
assert decode_expr(encode_expr(h)) == h
print(s.diff(h, z))
print(hpl_at_one((0, 0, 1)))        # zeta(3)
print(hpl_alternating_at_one((-1,)))  # log(2)
```

HPL expanded letters are `0, 1, -1`, with kernels `1/t`, `1/(1-t)`,
`1/(1+t)`. Thus `H_1(z)=-log(1-z)`. Constructors preserve symbolic words;
reductions are separate operations with documented domains.

## High-precision numerical values

```python
from hypra import evaluate_hpl

real = evaluate_hpl((1, 0), s.Rational(1, 3), dps=60)
complex_value = evaluate_hpl(
    (0, 1, -1), s.Rational(1, 4) + s.I/5,
    domain="principal-disk", dps=60,
)
print(real.value, real.tail_bound, real.stability_gap)
```

The default numerical domain is the real interval. Complex mode must be
selected explicitly and excludes its branch cut and unit-circle boundary.
Tail bounds concern truncation, not certified floating-point roundoff.
See [numerical conventions](hpl-numerics.md) and the separate
[radial-continuation API](hpl-continuation.md).

## Validate locally

```sh
uv run pytest -q
uv run ruff check .
uv run python scripts/validate.py --numeric --output-dir reports/local
```

The optional real-backend tests require an independently installed, licensed
Wolfram environment plus external HypExp/HPL. Skipped real-backend tests must
not be reported as a validated Mathematica comparison.
