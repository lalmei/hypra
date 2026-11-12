# Real classical weight-three reduction over {0,-1}

`classical_hpl(..., domain="0<z<1")` reduces all eight weight-three words
over `{0,-1}`, in addition to `{0,1}`. Words involving both `1` and `-1`
at weight three remain lossless HPL objects.

The shuffle identity depends only on words, so relabeling the alphabet in
its leading-letter/trailing-zero extraction gives a reduction to `00(-1)`
and `0(-1)(-1)`. This is an abstract relabeling of shuffle algebra, **not**
a substitution `z→-z` into branch-sensitive logarithmic formulas.
The first irreducible word is `-Li₃(-z)`. Write `L=log(1+z)` and
`u=1/(1+z)` for the second:

\[
 H_{0,-1,-1}(z)=\zeta(3)-Li_3(u)-L Li_2(u)
 +\tfrac12 L^2\log z-\tfrac13 L^3.
\]

Every logarithm and polylogarithm in this formula is real on `(0,1)`, with
`1/2<u<1`. No implicit cancellation of complex boundary values is needed.

To derive the formula, use the defining integral
`H_(0,-1,-1)(z)=½∫₀ᶻ log²(1+t) dt/t`. Substitution `v=1/(1+t)` gives

\[
 \frac12\int_u^1\log^2v\left(\frac1v+\frac1{1-v}\right)dv.
\]

A primitive for the second kernel is
`-log²v log(1-v)-2 log(v)Li₂(v)+2Li₃(v)`, verified by differentiation.
Evaluating its convergent upper limit and using `log u=-L` and
`log(1-u)=log z-L` yields the displayed identity. Its limit at zero is zero,
fixing the integration constant.

Tests check all eight exact differentiated recursions, the new formula's
origin limit and real polylogarithm arguments, and independent HPL numerical
values at `1/3` and `3/4` to 29 digits. The proof uses integral substitution,
endpoint constants, and shuffle identities, not fitted numerical values.
