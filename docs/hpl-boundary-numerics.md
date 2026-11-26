# Numerical real-cut boundary values

`evaluate_hpl_boundary` computes a declared upper or lower real-axis boundary
value using exact symbolic maps or Taylor transport to the exact real target.
No small imaginary displacement is substituted for the boundary limit.

```python
from hypra.hpl.boundary_numeric import evaluate_hpl_boundary

result = evaluate_hpl_boundary((0, 1), 2,
    region="x>1", boundary="upper", dps=60)
print(result.value)  # pi**2/4 + I*pi*log(2), numerically
print(result.stability_gap)
```

Both `region` and `boundary` are mandatory. Regions are `0<x<1`, `-1<x<0`,
`x>1`, `x<-1`; boundary is `upper` for `+i0` or `lower` for `-i0`. These are
principal boundary limits with **no extra winding**, not arbitrary paths.
The exact finite real input must satisfy the declared region; integers and
exact SymPy numbers are accepted, floats, strings, unknown inequalities and
symbolic free parameters are rejected. The singular endpoints `-1,0,1` are
excluded even for words with individually finite limits there.
For a finite value at exactly `-1`, use
[`hpl_at_minus_one`](hpl-alternating-endpoints.md) with the desired side.

## Method and branch checks

The existing [real-region continuation](hpl-real-continuation.md) transforms the
word algebraically to HPLs with arguments strictly between zero and one,
plus explicit powers of `I*pi` and endpoint constants. Those interior HPLs
use the convergent series evaluator. EulerSum/MZV constants use their
strict-weight Mobius recursion. No integration limit is interchanged by this
numerical composition and no tiny imaginary regulator is introduced.

Higher exterior weights use the triangular suffix differential equations.
The path starts at `1/4`, moves to `1/4 ± i/2`, then to `x ± i/2`, and
ends exactly at `x`. It stays on the selected side of every kernel pole.
Each Taylor step lies within a disk avoiding `0,±1`. This transport also
handles interior targets whose direct series exceeds its term budget.

For example, `H_1(2+i0)=+i*pi`, `H_0(-2+i0)=log(2)+i*pi`, and
`H_-1(-2+i0)=+i*pi`. Lower boundary values are their complex conjugates for
these real-kernel words. Tests independently check those logarithmic signs,
`Li2(2+i0)=pi²/4+i*pi*log(2)`, mixed-word shuffle identities and conjugacy.

## Accuracy diagnostics and limits

`HPLBoundaryNumericResult` records the complex value, region, boundary, path
convention, working/requested precision, observed precision-stability gap,
and `component_tail_bounds`. For symbolic-map evaluation that field contains **individual interior
HPL truncation bounds only**: they are not propagated through transformation
coefficients or cancellation, and are not a final output error bound.
Elementary calls can have zero mathematical truncation without zero roundoff.

The complete transformed expression is evaluated at `dps+20` and `dps+40`.
Taylor transport uses `dps+25` and `dps+45`, with a propagated truncation
estimate in `tail_bound` and an empty component-bound tuple. Failure of the
stability or resource checks raises an error. These checks do **not** certify floating-point
roundoff or all output digits using interval arithmetic.

All four regions support words through weight 24. Requested
precision is 2–800 decimal digits; component `max_terms` defaults to 100000.
Near singular endpoints the interior arguments can approach one and require
many terms. Symbolic map growth and constant recursion can also be costly.
No automatic dispatch from the ordinary real numerical API is introduced.
