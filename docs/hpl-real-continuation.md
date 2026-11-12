# Explicit real-axis boundary continuation

`analytic_continue_hpl(word, x, region=..., boundary=...)` covers four open
real regions: `"0<x<1"`, `"-1<x<0"`, `"x>1"`, and `"x<-1"`. The exact
strings are domain contracts supplied by the caller, not conditions inferred
from an otherwise unrestricted symbol. The mandatory boundary is `"upper"`
for `H_word(x+i0)` or `"lower"` for `H_word(x-i0)`.

The unit interval returns the unchanged HPL. The other regions produce exact
expressions in HPLs with arguments inside `(0,1)`, elementary constants, and
possibly unreduced alternating endpoint constants. The points `0,+1,-1` are
not included. This API describes standard boundary values, not arbitrary
winding paths or repeated monodromy around singular points.

## Composition

* For `-1<x<0`, set `t=-x`. Negation at the upper negative cut gives explicit
  `+i*pi` logarithms; reversing explicit `i` gives the lower boundary.
* For `x>1`, apply the inversion formula at `t=1/x` with the requested
  boundary unchanged.
* For `x<-1`, first write the negation identity at `x=-y`. Continue the
  resulting HPLs to `y>1`, then apply inversion at `t=1/y=-1/x`.
  Crucially, `x=-y` reverses the side: an upper boundary for `x` requires a
  **lower** boundary for `y`, and conversely. The explicit logarithmic
  constants from negation keep the sign selected for `x`.

Intermediate logarithms are converted to HPLs before mapping. Unknown
variable-dependent functions or unmapped variable coefficients are rejected.
Only real endpoint constants and explicit `i` appear in these compositions;
no unspecified conjugation assumptions are introduced.

Weight is bounded by 24 on the near intervals and six where inversion is
needed. Component algorithms may impose further resource limits. The caller
must not combine a region assertion with contradictory symbol assumptions.

Tests independently check every weight-one function in all four regions on
both boundary sides. In particular, `H_1(x)` has **positive** imaginary part
on the upper `x>1` boundary, whereas `H_0(x)` and `H_-1(x)` acquire positive
imaginary parts on their respective upper negative cuts. Further checks
verify the real dilogarithm at a negative argument beyond `-1` using its
classical inversion identity and exercise invalid region/boundary limits.
