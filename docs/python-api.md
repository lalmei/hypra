# Python API reference

The package root `hypra` exports the public functions and classes below.
The custom backend protocol lives in `hypra.api`. The signatures and
Google-style docstrings are rendered from the implementation by
`mkdocstrings-python`. See [API examples](api-examples.md) for runnable
examples and [capabilities and limits](scope.md) for supported domains.

An `IntegerWord` is a tuple, list, or SymPy tuple of Python or SymPy integers.
An `ExactScalar` is a SymPy expression or a Python integer. Individual functions
describe any extra requirements on letters, exactness, and argument domains.

## Expansion and backends

::: hypra.api.expand_epsilon

::: hypra.api.backend_info

::: hypra.api.ExpansionResult
    options:
      members: [expression, coefficient, latex, to_dict]

::: hypra.api.ExpansionBackend

::: hypra.native.NativeConfig

::: hypra.expansion_library.ExpansionLibrary

::: hypra.native.NativeBackend

::: hypra.config.BackendConfig
    options:
      members: [from_env, resolved_executable, resolved_cache_dir]

::: hypra.hypergeometric

Use `hypergeometric` when every supplied parameter needs to remain visible
for validation. See the [exceptional-parameter rules](scope.md#preserving-exceptional-input-parameters).

`NativeConfig.argument_domain` selects the [interior branch](scope.md#principal-complex-arguments).
`NativeConfig.unit_endpoint` selects the [unit-endpoint policy](endpoints.md).
`NativeConfig.unit_excess_shift` selects a lower parameter to shift at a
[fixed excess pole](endpoints.md#fixed-excess-poles).
See the [runnable examples](api-examples.md) for these settings.

## Exact transport

::: hypra.serialization.encode_expr

::: hypra.serialization.decode_expr

::: hypra.serialization.SerializationError

## Expansion errors

::: hypra.errors.HypraError

::: hypra.errors.UnsupportedExpressionError

::: hypra.errors.IncompleteExpansionError

::: hypra.backend.BackendError

::: hypra.backend.BackendUnavailableError

::: hypra.backend.BackendTimeoutError

::: hypra.backend.BackendProtocolError

## Symbolic functions and notation

::: hypra.symbolic.HPL
    options:
      members: [__new__, from_compressed, word, argument]

::: hypra.hpl.notation.expand_indices

::: hypra.hpl.notation.compress_indices

::: hypra.hpl.notation.pm_to_letters

::: hypra.hpl.notation.letters_to_pm

::: hypra.hpl.notation.hpl_pm

::: hypra.nielsen.NielsenPolylog
    options:
      members: [__new__, n, p, argument, as_hpl]

::: hypra.nielsen.hpl_to_nielsen

::: hypra.nielsen.nielsen_to_hpl

::: hypra.nielsen.reduce_nielsen

## Algebra and integration

::: hypra.hpl.algebra.shuffle_expand

::: hypra.hpl.algebra.primitive_hpl

::: hypra.hpl.algebra.integrate_from_zero

::: hypra.hpl.algebra.hpl_series

::: hypra.hpl.classical.classical_hpl

::: hypra.hpl.reduction.extract_hpl_logs

::: hypra.hpl.reduction.lyndon_factors

::: hypra.hpl.reduction.lyndon_reduce

::: hypra.hpl.taylor.hpl_taylor

## Endpoints and constants

::: hypra.symbolic.MZV
    options:
      members: [__new__, indices]

::: hypra.symbolic.EulerSum
    options:
      members: [__new__, indices, colors]

::: hypra.hpl.endpoints.hpl_at_one

::: hypra.hpl.endpoints.reduce_mzv

::: hypra.hpl.alternating.hpl_alternating_at_one

::: hypra.hpl.alternating.reduce_euler_sum

::: hypra.hpl.regularization.hpl_endpoint_polynomial

::: hypra.hpl.regularization.shuffle_regularized_hpl_at_one

::: hypra.hpl.regularization.hpl_expression_endpoint_polynomial

::: hypra.hpl.stuffle.stuffle_expand_constants

::: hypra.hpl.special.hpl_at_i

::: hypra.hpl.series_one.hpl_series_at_one

::: hypra.hpl.limits.hpl_limit_at_one

::: hypra.hpl.named_constants.HPLs6

::: hypra.hpl.named_constants.HPLs7a

::: hypra.hpl.named_constants.HPLs7b

::: hypra.hpl.named_constants.HPLs8a

::: hypra.hpl.named_constants.HPLs8b

::: hypra.hpl.named_constants.HPLs8c

::: hypra.hpl.named_constants.HPLs8d

## Argument maps

::: hypra.hpl.arguments.transform_hpl_arguments

::: hypra.hpl.transforms.negate_hpl

::: hypra.hpl.transforms.reflect_hpl

::: hypra.hpl.maps.square_hpl

::: hypra.hpl.maps.fractional_hpl

::: hypra.hpl.mobius.mobius_hpl

::: hypra.hpl.inversion.invert_hpl

::: hypra.hpl.real_continuation.analytic_continue_hpl

## Numerical evaluation

::: hypra.hpl.numeric.evaluate_hpl

::: hypra.hpl.numeric.hpl_evalf

::: hypra.hpl.numeric.HPLNumericResult

::: hypra.hpl.numeric.HPLNumericalConvergenceError

::: hypra.hpl.radial_continuation.continue_hpl

::: hypra.hpl.radial_continuation.HPLContinuationResult

::: hypra.hpl.boundary_numeric.evaluate_hpl_boundary

::: hypra.hpl.boundary_numeric.HPLBoundaryNumericResult

::: hypra.hpl.constants.evaluate_mzv

::: hypra.hpl.constants.evaluate_euler_sum

::: hypra.hpl.constants.ConstantNumericResult

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

::: hypra.hypexp_integrals.beta_log_moment

::: hypra.hypexp_integrals.hypergeometric_log_moment
