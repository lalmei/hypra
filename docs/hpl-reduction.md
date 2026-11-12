# Exact HPL logarithm extraction

`extract_hpl_logs(expr)` removes leading 1 letters and trailing 0 letters by
shuffle identities, retaining explicit powers of `-log(1-z)` and `log(z)`.
For example H(1,0;z) = -log(z)log(1-z) - Li₂(z).

This preserves the argument and does not choose a new continuation path.
Repeated identical letters reduce to the corresponding weight-one HPL power
divided by the factorial. General remaining words stay exact; no claim of a
minimal transcendental basis is made. Each shuffle retains multiplicity.

The recursive trailing-zero reduction decreases the number of trailing zeros
at fixed weight; leading-one removal similarly decreases leading ones. Lower
weight factors terminate the recursion. Identities are independently tested
using local exact HPL series, including mixed-sign words.

## Table-free shuffle basis

`lyndon_reduce(expr)` expresses words through weight eight as polynomials in
Lyndon-word HPLs, ordered by the alphabet −1 < 0 < 1. Duval's algorithm gives
the unique nonincreasing Chen–Fox–Lyndon factorization of each word. Shuffling
those factors produces the original word as its lexicographically greatest
term, with a nonzero integer multiplicity. Solving this identity recursively
therefore decreases lexicographic order at fixed weight and terminates.
The implementation checks triangularity instead of assuming it silently.

Repeated factors retain their shuffle multiplicities, for example
H(1,1,1;z)=H(1;z)³/6. This is a basis of the formal shuffle algebra, not a
claim that its evaluations are algebraically independent after special values,
argument transformations, or additional functional relations. No upstream
reduction tables are copied. Word/weight limits fail explicitly.

## Local primitives versus definite integrals

`primitive_hpl(expr,z)` returns a local antiderivative up to an arbitrary
constant, including cases singular at zero. Its rational coefficients may
have poles at 0,±1; products are shuffle-linearized. `integrate_from_zero`
uses the same primitive but additionally proves that the combined origin
limit is finite and subtracts it. A primitive for `1/z` is therefore valid,
while its integral from zero is rejected. Neither API invents a regularized
endpoint value or a new contour prescription.
