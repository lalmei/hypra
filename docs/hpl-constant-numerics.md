# Numerical Euler sums and multiple zeta values

`evaluate_euler_sum(indices, colors, dps=50)` and
`evaluate_mzv(indices, dps=50)` evaluate convergent endpoint constants through
weight six. They return a `ConstantNumericResult` with `.value`,
`.stability_gap`, requested/working precision, method and sum convention.
Indices and colors have the exact conventions of [alternating endpoints](hpl-alternating-endpoints.md).
Divergent constants are rejected rather than regularized.

```python
from hypra.hpl_constants import evaluate_euler_sum, evaluate_mzv

print(evaluate_euler_sum((1, 2), (-1, -1), dps=60).value)
print(evaluate_mzv((4, 2), dps=60).value)
```

## Why endpoint summation is unnecessary

Given a colored sum, let `a_j` be the cumulative product of its first `j`
colors. Build its word from blocks `0**(s_j-1), a_j`, and let
`f=prod(a_j)`. The endpoint identity is `H_word(1)=f*EulerSum`.
A convergent sum gives a word beginning in `0` or `-1`.

The exact [Mobius transformation](hpl-mobius.md) has the form

\[
H_w\!\left(\frac{1-z}{1+z}\right)=C_w+R_w(z),\qquad C_w=H_w(1).
\]

Its differential recursion integrates a transformed suffix of weight one
less. Therefore every endpoint constant in `R_w` has strictly smaller weight.
Subtracting the identical exact `C_w` expression before numerical evaluation
cancels current-weight constants algebraically. At `z=1/2` this gives

\[
C_w=H_w(1/3)-R_w(1/2).
\]

The evaluator recursively replaces lower-weight EulerSum/MZV atoms and uses
the convergent HPL numerical series at `1/3` and `1/2`. It explicitly checks
that no surviving constant has weight greater than or equal to the current
weight. An unsuccessful reduction fails rather than entering a circular
numeric definition. Exact elementary/Newton/MZV reductions are used first.
All evaluation stays on the real interval; no numerical analytic
continuation or branch selection is required.

## Accuracy and scope

Two calculations at `dps+20` and `dps+40` must agree more closely than the
requested relative precision. HPL component calls impose their stricter
analytic series-tail targets. The reported stability gap is evidence, **not
an interval-certified error bound for the final constant**: cancellation,
roundoff and special-function evaluation are not rigorously interval-enclosed.
A term-budget or precision-stability failure raises an error, with no partial
success result. The evaluator preserves global mpmath precision.

Supported weights are at most six, requested precision 2–900 decimal digits,
and HPL `max_terms` defaults to 100000. Exact Möbius expressions can grow
substantially with weight, so weight-six calls may be slower than elementary
reductions. Explicit `MZV(...).evalf(50)` and `EulerSum(...).evalf(50)` use this
engine through weight six, with guarded digits. Construction remains
symbolic; constants above the supported weight remain unevaluated. Tests include independent one-dimensional quadrature and
convergent finite-cutoff stuffle identities, as well as the symbolic strict
weight-decrease invariant.
