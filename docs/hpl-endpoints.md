# HPL endpoints and exact multiple zeta values

`hpl_at_one(word)` evaluates a discrete HPL word at **z = 1**, taking the
limit from the real interval `0 < z < 1`. This first endpoint implementation
supports the alphabet `{0, 1}` only. It is not a generic limit engine for sums
with singular rational prefactors; cancellations between divergent HPLs must
not be inferred by evaluating their endpoints independently.

```python
from hypra import MZV, hpl_at_one, reduce_mzv

hpl_at_one((0, 1))                 # pi**2/6
hpl_at_one((1, 0))                 # -pi**2/6 (finite, despite leading 1)
hpl_at_one((0, 1, 0))              # -2*zeta(3)
hpl_at_one((0, 0, 1, 0, 1), reduced=False)  # MZV((3, 2))
reduce_mzv(MZV((2, 1)))            # zeta(3)
```

## Conventions and convergence

The exact symbolic class uses strict descending summation:

\[
\zeta(s_1,\ldots,s_d)
=\sum_{n_1>\cdots>n_d\ge1}\frac{1}{n_1^{s_1}\cdots n_d^{s_d}},
\qquad s_1\ge2,\quad s_j\in\mathbb Z_{>0}.
\]

Divergent indices such as `(1,2)` are rejected; they are not assigned a
regularized value. Empty indices equal one. Other MZVs remain symbolic unless
an explicitly proved identity below applies. MZV objects have exact JSON AST
round-tripping and LaTeX output; no numerical fit supplies their values.

An admissible HPL word begins with `0` and ends with `1`. Grouping it into
runs `0^(s1-1) 1 ... 0^(sd-1) 1` gives the corresponding MZV indices. Nonempty
all-zero words have endpoint zero, and the empty HPL word has endpoint one.

Trailing zeros are removed using the shuffle product with `H_0(z)=log(z)`.
At `z=1`, `log(z) H_w(z)` tends to zero: HPLs on this alphabet grow at most as
powers of `log(1-z)`. Solving the shuffle identity removes one trailing zero
without dropping a divergent summand. For example,

\[
H_{1,0}(z)=\log z\,H_1(z)-H_{0,1}(z),\qquad
H_{1,0}(1)=-\zeta(2).
\]

Words starting with `1` and containing a second `1` have a nonzero or
divergent tail at the endpoint, and are rejected as divergent. A leading `1`
followed only by one or more zeros is instead integrable. The implementation
never assigns a finite value to a divergent endpoint.

## Proved reductions, not a table lookup

MZV duality reverses an admissible HPL word and exchanges `0` with `1`.
It follows by substituting `t -> 1-t` and reversing the iterated-integration
simplex. A shorter dual representative can reduce directly to depth one;
for example `zeta(2,1)=zeta(3)` and `zeta(2,1,1)=zeta(4)`.

The weight-four identities are obtained by comparing shuffle and stuffle:

\[
\zeta(2)^2=2\zeta(2,2)+\zeta(4)
=2\zeta(2,2)+4\zeta(3,1).
\]

Together with `zeta(2)^2 = 5 zeta(4)/2`, these give

\[
\zeta(3,1)=\tfrac14\zeta(4),\qquad
\zeta(2,2)=\tfrac34\zeta(4).
\]

No general MZV basis reduction, alternating endpoint sums, or transcendence
claim is implemented. Endpoint/duality combinatorics are limited to weight 24;
higher-weight MZVs remain exact symbolic objects. Depth-one reduction is not
subject to that combinatorial limit. The optional Wolfram comparison bridge
does not transport MZV heads; this is a native-port feature.
