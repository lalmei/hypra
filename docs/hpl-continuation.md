# Explicit numerical radial HPL continuation

```python
import sympy as s
from hypra.hpl_continue import continue_hpl

answer = continue_hpl((0, 1), 3*s.I, dps=40,
                      domain="principal-cut-plane", path="radial")
print(answer.value)          # numerical Li_2(3i)
print(answer.tail_bound)     # absolute mathematical truncation estimate
print(answer.stability_gap)  # independent working-precision comparison
```

Both `domain` and `path` must be specified. The only implemented path is the
straight radial segment from a point in the principal origin disk to the
requested target. Real points are allowed only in (0,1); zero, ±1 and real
cuts (-infinity,0] and [1,infinity) are rejected. Nonreal targets may lie outside
the unit disk. No endpoint boundary value, winding path, or arbitrary contour
is inferred. This is numerical continuation, not exact symbolic continuation.
The resulting branch follows the explicitly selected path from the principal
origin germ; it is not a branch-independent value.

## Algorithm and truncation estimate

For word w, collect its suffix functions and the empty word H_empty=1 in a
state vector y. The defining HPL derivatives give a triangular system

\[
 y'_j=f_{w_j}(z)y_{j+1},\quad
 f_0=1/z,\quad f_1=1/(1-z),\quad f_{-1}=1/(1+z).
\]

Start at c=z min(1,1/(4|z|)). Each suffix is initialized with the existing
`evaluate_hpl(..., domain="principal-disk")`, including its mathematical
series tail estimate. A local step h obeys |h|<=d/4, where
`d=min(|c|,|c-1|,|c+1|)`. Taylor coefficients are obtained by exact-form kernel
power-series convolution, evaluated at high floating-point precision.

On the complex disk of radius R=d/2, the matrix infinity norm is at most
1/(d-R). Thus the exact ODE solution starting from the current approximate
state has norm bounded by

\[
 M=\|\widetilde y(c)\|_\infty\exp(R/(d-R)).
\]

Cauchy's coefficient estimate bounds the degree-N Taylor tail at h by

\[
 M\frac{(|h|/R)^{N+1}}{1-|h|/R}.
\]

The accumulated initial/truncation error E is propagated separately, never
omitted or absorbed into the approximate state's norm:

\[
 E_{new}\le E\exp(|h|/(d-|h|))
       +M\frac{(|h|/R)^{N+1}}{1-|h|/R}.
\]

This follows from the same matrix norm bound along the step segment. The
maximum initial suffix-tail bound seeds E. The step order is selected from
the local estimate, and the final propagated estimate must meet the requested
scaled absolute tolerance `10**(-dps)*max(1,abs(value))`; otherwise the call
raises. It does not silently return an unchecked partial continuation.

## What the estimate does not certify

These are mathematical truncation estimates evaluated with floating-point
arithmetic, not outward-rounded interval bounds. They do not include all
roundoff errors, the finite decimal conversion of the initial point, or
uncertainty already present in a float input. A second
run at higher working precision must agree to the requested tolerance, but
that stability check is evidence, not a rigorous certification of digits.
Prefer exact SymPy numerical arguments. The returned bound is absolute;
near a zero it is not a guarantee of relative significant digits.

The implementation has explicit limits on precision, step count, local
Taylor degree, origin-series terms and HPL weight. Targets very close to a
real cut can require too many steps and fail. It does not change the exact
native expansion API's accepted branch domain.

Focused tests compare logs and classical polylogarithms beyond the unit disk,
a mixed trailing-log shuffle identity, lower-half-plane values and rejection
of ambiguous domains/endpoints. Existing exact HPL conventions remain unchanged.
