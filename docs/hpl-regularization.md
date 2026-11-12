# Explicit shuffle-regularized endpoints

These APIs construct an endpoint polynomial or its constant term. They do
**not** declare a divergent HPL to have a finite limit.

```python
import sympy as s
from hypra.hpl_regularized import (
    hpl_endpoint_polynomial, shuffle_regularized_hpl_at_one,
)
T = s.Symbol("T")
hpl_endpoint_polynomial((1, 0, 1), T)  # zeta(2)*T - 2*zeta(3)
shuffle_regularized_hpl_at_one((1,))   # 0, a regularized constant, NOT a limit
```

## Fixed convention

The approach is real `z -> 1-`, and `T=-log(1-z)=H_1(z)`. For an expanded
word over `0,1,-1`, return the exact polynomial `P_w(T)` such that
`H_w(z)-P_w(-log(1-z))` tends to zero. The discarded remainder also vanishes
when multiplied by any fixed power of `T`. The shuffle-regularized constant
is **by definition** `P_w(0)`.

There is no complex approach direction, implicit branch choice, subtraction
scale, or stuffle regularization in this API. The existing convergent
endpoint functions still reject genuine divergences. `expand_epsilon` does
not call these regularization helpers implicitly.

## Derivation and algebra

If a word does not begin with `1`, use its convergent alternating endpoint.
For leading ones, the exact shuffle identity is

\[
H_1(z)H_v(z)=\sum_w m_w H_w(z).
\]

Isolate the original word `w=(1,v)` with its shuffle multiplicity. Every
other term has fewer leading ones. Replace `H_1` by `T`, recursively replace
the remaining HPLs by their endpoint polynomials, and divide by that exact
multiplicity. All-one words give `T^n/n!`; the empty word gives one.

Regular-singular endpoint expansions have only nonnegative powers of `1-z`
times finite logarithmic powers. After the logarithmic polynomial is
removed, multiplication by any fixed power of `T` still tends to zero. Thus
the shuffle step preserves the polynomial asymptote, and products obey
`P_u P_v = sum_w m_w P_w`. Taking constant terms therefore respects shuffle
products, but does not imply an unrelated stuffle regularization convention.

This also handles convergent words that happen to begin with `1`, for example
`H(1,0;1)=-zeta(2)`: the apparent logarithmic coefficient vanishes exactly.

The implementation uses existing shuffle multiplicities and convergent
endpoint constants, not copied special-value tables. Its explicit word
limit is 12; general alternating constants remain symbolic. This is not an
endpoint expansion to positive powers of `1-z`, nor a general limit engine
for arbitrary rational-HPL expressions.

## References

The kernel and logarithm-extraction conventions follow D. Maître,
[HPL, a Mathematica implementation of the harmonic polylogarithms](https://arxiv.org/abs/hep-ph/0507152),
and the project's [source and algorithm map](upstream-source.md).
The constant-term convention above is specified explicitly rather than
claiming equivalence to every upstream regularization or analytic-continuation
option.
