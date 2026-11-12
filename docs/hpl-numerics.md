# High-precision native HPL evaluation

`hpl_evalf(word, argument, dps=50)` evaluates a word in `{0,1,-1}` on the
**real interval `0 < argument < 1`**, using the principal real HPL branch.
It returns a high-precision `mpmath` real. This is an approximate numerical
API, separate from the package's exact symbolic constructors and AST.

```python
from sympy import Rational
from hypra import hpl_evalf, evaluate_hpl

value = hpl_evalf((0, 1), Rational(1, 3), dps=80)
check = evaluate_hpl((1, 0, -1), "0.4", dps=60)
print(check.value)
print(check.terms, check.tail_bound, check.stability_gap)
```

`evaluate_hpl` returns the value plus its method, requested and working
precision, number of series terms, analytic truncation bound and observed
precision-stability gap. Repeated-letter words use logarithms; words of zeros
followed by `1` or `-1` use classical polylogarithms. Other words, including
mixed trailing-zero logarithmic words, use the convergent coefficient
recurrence below. There is no nested numerical quadrature.

## What the accuracy checks do and do not establish

For the series method, the omitted mathematical tail is bounded using a
majorant of **all subsequent coefficients**, not merely a small last term.
That bound is evaluated at guarded numerical precision. The complete
calculation is repeated at `dps+20` and `dps+40` working decimal digits, with
a stricter truncation target than the requested precision. A failure to reach
the target within `max_terms`, or disagreement between the two working
precisions, raises `HPLNumericalConvergenceError`; no partial value is returned.

**This is not interval certification of floating-point roundoff.**
`tail_bound` concerns the series truncation only. A zero tail bound for a
logarithm/classical-polylogarithm call does not claim that mpmath incurred no
roundoff. `stability_gap` is evidence of precision stability, not a proof of
all output digits. Tests compare independent elementary/polylog identities,
single-integral formulas and shuffle identities.

For an intended decimal or rational point, pass a decimal string or a SymPy
`Rational`. Python floats and finite-precision mpmath arguments are evaluated
at their stored numerical value; increasing `dps` does not recover missing
input digits. The evaluator does not change the caller's global mpmath
precision. Numeric strings are parsed only as decimal numbers, not as code.

## Coefficient recurrence and tail bound

Write the real germ as

\[
H_w(z)=\sum_{n\ge0}z^nP_{w,n}(\log z).
\]

For a trailing block of `d` zeros, `P` has degree at most `d` in its logarithm.
Pure-zero suffixes have the exact base `log(z)^d/d!`. Integration solves

\[
(n+\partial_L)Q_n(L)=P_n(L).
\]

A leading zero applies that integration to the same tail coefficient. A
leading `1` or `-1` first accumulates the preceding tail coefficients with
kernel signs `1` or alternating signs, respectively. Thus the work is linear
in the number of series terms, word weight and log degree, with memory
independent of the term count.

For `R=max(1,|log z|)` and

\[
C_d=\sum_{j=0}^{d}\frac{d!}{(d-j)!},\qquad
B=C_d^{|w|-d}\frac{R^d}{d!},
\]

the coefficient-polynomial norm is uniformly bounded by `B`. Indeed, the
integration operator has norm at most `C_d/n`, while a nonzero kernel sums
at most `n` previous coefficients. Therefore the omitted tail after term `N`
satisfies the mathematical bound

\[
\left|\sum_{n>N}z^nP_{w,n}(\log z)\right|
\le B\frac{z^{N+1}}{1-z}.
\]

The implementation includes extra numerical slack when computing this
majorant; it still makes no claim of rigorous interval-rounded arithmetic.

## Deliberate limits

- Weight at most 24; requested `dps` from 2 to 1000.
- Default `max_terms=100000`; configurable up to one million.
- General words near one may need many terms and can exhaust that budget.
- Endpoints, complex arguments and other branches are rejected. Use
  `hpl_at_one` for the supported exact endpoint cases.
- There is no automatic `.evalf()` hook on the symbolic HPL class yet.
- Complex unit-disk continuation is a separate future increment, not silently
  enabled by this API.
