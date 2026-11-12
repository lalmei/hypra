# Alternating endpoint constants

`hpl_alternating_at_one(word, reduced=True)` evaluates the finite limit of an
expanded HPL word over `{0,1,-1}` as its real argument increases to one.
It accepts weight up to 24, including mixed trailing zeros. It does not
assign finite parts or regularized values to divergent words.

The exact symbolic `EulerSum(indices, colors)` uses the explicit convention

\[
\operatorname{ES}(s_1,\ldots,s_d;\sigma_1,\ldots,\sigma_d)
=\lim_{N\to\infty}\sum_{N\ge n_1>\cdots>n_d\ge1}
\prod_{j=1}^d\frac{\sigma_j^{n_j}}{n_j^{s_j}},\qquad \sigma_j\in\{-1,1\}.
\]

Indices are positive integers. Convergence requires `s1 > 1` or first color
`-1`. In the latter case the sum may be conditionally convergent: the
outer-index cutoff above is part of the definition, not permission to reorder
terms. In particular `ES((1,),(-1,)) = -log(2)`. Empty indices/colors give one.
The constructor preserves the constant symbolically; AST transport records
both exact lists, and LaTeX shows indices and colors in separate rows.
The optional Wolfram bridge does not currently decode this native AST head.

For a word ending in a nonzero letter, let its successive nonzero letters be
`a1,...,ad`, and let `sj-1` count the zeros immediately before `aj`. Then

\[
H_w(1)=(a_1\cdots a_d)\operatorname{ES}
(s_1,\ldots,s_d; a_1,a_1a_2,\ldots,a_{d-1}a_d).
\]

This follows by expanding the kernels `1/(1-t)` and `1/(1+t)` into geometric
series on the open unit interval and taking the convergent endpoint (Abel
limit in conditional cases). Trailing zeros are removed using the shuffle
product with `H_0(z)=log(z)`. Its product with any remaining word vanishes at
one because endpoint growth is at most a finite power of `log(1-z)`.

A leading `1` is divergent unless its entire tail consists of zeros. All HPL
kernels are positive on the real interval; a tail containing a nonzero letter
has a fixed sign and nonzero endpoint magnitude, so it cannot cancel the
leading `1/(1-z)` singularity. Empty words give one and nonempty all-zero
words give zero. This test rejects genuine divergences without pretending
that every leading `1` word diverges: `H(1,0;1) = -zeta(2)` is finite.

Opt-in reductions are deliberately small and derived, not copied tables:

- All-positive colors become ordinary MZVs, using existing proved reductions.
- Depth-one negative colors use the Dirichlet eta identity
  `ES(n;-1)=-(1-2**(1-n))*zeta(n)` for `n>1`, and `-log(2)` for `n=1`.
- The repeated `-1` HPL word gives `log(2)**d/d!` by repeated integration.

Other constants remain exact EulerSum objects. `reduced=False` retains even
these elementary constants in their lossless sum representation. This is
not a complete basis reduction for alternating multiple zeta values.

## Repeated indices and colors: finite-cutoff Newton identities

For `ES((s,)*d,(c,)*d)`, let `x_n=c**n/n**s` and let `e_d(N)` be the
coefficient of `t**d` in the finite polynomial `prod(n=1..N)(1+t*x_n)`.
Its coefficient is exactly the finite strict nested sum. The logarithmic
derivative of this polynomial gives Newton's identity

\[
e_0(N)=1,\qquad
d e_d(N)=\sum_{k=1}^{d}(-1)^{k-1}e_{d-k}(N)p_k(N),\qquad
p_k(N)=\sum_{n=1}^{N}\frac{(c^k)^n}{n^{ks}}.
\]

Every power sum converges for a valid repeated input: either `s>1`, or
`s=1,c=-1` with a convergent alternating first power and absolutely
convergent higher powers. Taking limits in this finite algebraic identity
therefore yields an exact recursive reduction to depth-one constants.
There is no rearrangement of a conditionally convergent infinite series.
For example,

\[
\operatorname{ES}(1,1;-1,-1)=\frac{\log^2 2-\zeta(2)}{2}.
\]

The implementation applies this recurrence at arbitrary finite depth,
not a table of special values. Tests include the explicit depth-two and
-three formulas and independently accumulated finite nested sums.
