# Python API reference

These signatures and docstrings come from the actual Python implementation,
using `mkdocstrings-python`. For complete mathematical conditions, read
[scope](scope.md) and the linked algorithm chapters. Generated API descriptions
are not a claim of additional parameter or branch coverage.

::: hypra.api.expand_epsilon

::: hypra.api.backend_info

::: hypra.api.ExpansionResult

::: hypra.native.NativeConfig

::: hypra.hypergeometric

Use `hypergeometric` when every supplied parameter needs to remain visible
for validation. See the [exceptional-parameter rules](scope.md#preserving-exceptional-input-parameters).

`NativeConfig.argument_domain` selects the [interior branch](scope.md#principal-complex-arguments).
`NativeConfig.unit_endpoint` selects the [unit-endpoint policy](endpoints.md).
`NativeConfig.unit_excess_shift` selects a lower parameter to shift at a
[fixed excess pole](endpoints.md#fixed-excess-poles).
See the [runnable examples](api-examples.md) for these settings.

::: hypra.symbolic.HPL

::: hypra.hpl_classical.classical_hpl

::: hypra.hpl_regularized.hpl_endpoint_polynomial

::: hypra.hpl_regularized.shuffle_regularized_hpl_at_one

::: hypra.hpl_stuffle.stuffle_expand_constants

::: hypra.hpl_regularized.hpl_expression_endpoint_polynomial

::: hypra.hpl_special.hpl_at_i

::: hypra.hpl_taylor.hpl_taylor

::: hypra.hpl_boundary_numeric.evaluate_hpl_boundary

::: hypra.hpl_boundary_numeric.HPLBoundaryNumericResult

::: hypra.nielsen.NielsenPolylog

::: hypra.nielsen.hpl_to_nielsen

::: hypra.nielsen.nielsen_to_hpl

::: hypra.nielsen.reduce_nielsen

::: hypra.hpl_series_one.hpl_series_at_one

::: hypra.hpl_limits.hpl_limit_at_one


## Definite logarithmic integrals

`beta_log_moment(n, m, p)` integrates
\(u^p \log^n(u) \log^m(1-u)\) from zero to one. It normalizes the
rational-HPL primitive at zero and combines its finite endpoint limit at
one. The powers are explicit integers with `n, m >= 0`, `n + m <= 7`,
and `-32 <= p <= 64`. Convergence at zero requires `p + m >= 0`.

`hypergeometric_log_moment(p, n, m, k, q, z)` adds the factor
\(\log^k(1-zu)/(zu-1)^q\). All five powers are nonnegative integers.
It supports `z=0`, convergent cases at `z=1` (`n >= q`), and
real interior arguments `0<z<1`, including products of all three
logarithms. A symbolic interior argument is interpreted on that interval. A provably
negative real argument is reflected to `z/(z-1)` by substituting `u=1-v`.
The logarithmic weight `n+m+k+int(q>0)` is at most five, with
`p <= 64` and `q <= 32`. The interior mixed-log recurrence
differentiates the Euler hypergeometric equation in the three log
parameters and fixes the solution at the origin. Integration by parts
raises the polynomial degree, and differentiation raises the pole order.
Both APIs return exact SymPy expressions and raise
`UnsupportedExpressionError` for unsupported cases or divergent integrals.

```python
from sympy import Symbol
from hypra import beta_log_moment, hypergeometric_log_moment

z = Symbol("z")
beta_log_moment(0, 1, -1)                       # -pi**2/6
hypergeometric_log_moment(0, 0, 0, 0, 1, z)    # -HPL((1,), z)/z
```
