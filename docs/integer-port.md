# Integer-parameter shifts for the Python port

This note derives the contiguous shifts from the defining hypergeometric series.
It is an independent mathematical specification, not a literal translation of
upstream Mathematica source. The source inspection map and attribution are in
[upstream-source.md](upstream-source.md).

Let

\[
F(\mathbf a;\mathbf b;z)=\sum_{n\ge0} C_nz^n,\qquad
C_n=\frac{\prod_i(a_i)_n}{\prod_j(b_j)_n\,n!},\qquad
\theta=z\frac{d}{dz}.
\]

Initially take \(|z|<1\), no denominator parameter a nonpositive integer, and
parameters in a neighborhood where the expressions below exist. Identities then
continue meromorphically in parameters and analytically in z along a consistently
chosen continuation path. They do not independently authorize a branch change or
an endpoint integration at z=1.

## 1. Raise an upper parameter

The Pochhammer ratio is
\((a+1)_n/(a)_n=(a+n)/a\), giving

\[
F(a+1,\widehat{\mathbf a};\mathbf b;z)
 =\frac{\theta+a}{a}F(a,\widehat{\mathbf a};\mathbf b;z).
\]

For target \(a=A+\alpha\epsilon\), \(A\ge1\), start at
\(1+\alpha\epsilon\) and raise A-1 times. Every divisor has a positive integer
constant term, so these shifts consume **no additional epsilon order**. Using
\(\alpha\epsilon\) as the base instead introduces an artificial division by
\(\epsilon\); do not do that when the base with integer part 1 is available.
If \(a=0\) exactly, the displayed division is invalid: the value of F at a=0
alone does not determine its value at a=1.

## 2. Lower a lower parameter

Since \((c)_n/(c-1)_n=(c+n-1)/(c-1)\),

\[
F(\mathbf a;c-1,\widehat{\mathbf b};z)
 =\frac{\theta+c-1}{c-1}F(\mathbf a;c,\widehat{\mathbf b};z).
\]

Starting at \(c=1+\beta\epsilon\), the first step to
\(\beta\epsilon\) divides by \(\beta\epsilon\). All subsequent steps down to
an integer part -m have nonzero constant divisors. Thus each lower parameter
crossing integer part 0 costs at most **one additional epsilon coefficient**,
provided \(\beta\ne0\). This is an intermediate-computation budget, even if
numerator zeros remove the pole from the final result. If \(\beta=0\), the
target has an exact nonpositive lower parameter and requires a separate explicitly
defined terminating/singular convention; the generic routine must reject it.
For symbolic beta, nonzero must be assumed or proved, not guessed.

A concrete order trap uses the z coefficient of
\(F={}_2F_1(1,1;1+\epsilon;z)\), which is \(1/(1+\epsilon)\).
The shifted coefficient is exactly \(1/\epsilon\). But using F only through
\(\epsilon^3\) gives

\[
\frac{1+\epsilon}{\epsilon}
 (1-\epsilon+\epsilon^2-\epsilon^3)
 =\frac1\epsilon-\epsilon^3.
\]

The omitted \(\epsilon^4\) coefficient is needed to cancel the spurious
\(\epsilon^3\) term. Truncation order must be planned **before** shifting.

## 3. Raise a lower parameter

The coefficient ratio is \((c)_n/(c+1)_n=c/(c+n)\). For \(\Re c>0\),

\[
F(\mathbf a;c+1,\widehat{\mathbf b};z)
 =c z^{-c}\int_0^z t^{c-1}
      F(\mathbf a;c,\widehat{\mathbf b};t)\,dt.
\]

For positive target integer part C, start at \(1+\beta\epsilon\) and apply
this C-1 times. For epsilon in a sufficiently small disk, \(\Re c\) stays
bounded below by a positive number. The base series is analytic near t=0 and
has value 1 there. Its epsilon derivatives introduce finite powers of log(t),
which remain integrable against \(t^{\delta-1}\) for any fixed \(\delta>0\).
On compact subintervals of 0<z<1, these bounds justify termwise epsilon expansion.
Expanding \(t^{\beta\epsilon}\) and \(z^{-\beta\epsilon}\) requires retaining
all terms through the requested order but introduces **no extra pole order**.

The integral is an antiderivative with its lower-limit constant fixed by the
analytic hypergeometric germ F(0)=1, not an indefinite primitive with its constant
silently discarded. Apparent log(z) or negative powers from intermediate HPL
forms must cancel as required by that normalization. For complex z, the path,
logarithm and HPL branches must agree. This local argument says nothing about
interchanging a further integral over z up to the singular endpoint 1.

## 4. Negative upper integers: safe finite-prefix decomposition

Naively inverting the upper-raising operator gives a zero-endpoint integral that
diverges for nonpositive target integer parts. Do not use that primitive without
subtractions and analytic-continuation data. A safer exact alternative is the
finite-prefix/tail identity, for any positive integer N:

\[
F(\mathbf a;\mathbf b;z)
 =\sum_{n=0}^{N-1}C_nz^n
  +C_Nz^N\,{}_{p+1}F_{q+1}
   \left(1,\mathbf a+N;N+1,\mathbf b+N;z\right).
\]

Proof: set n=N+k in the tail and use
\((a)_{N+k}=(a)_N(a+N)_k\) and
\((N+k)!=N!(N+1)_k\); the new upper parameter 1 restores the standard 1/k!.
Choose N so the troublesome integer parts become positive. This splits the
finite terms that would otherwise collide with denominator shifts before any
epsilon expansion. The tail generally increases p and q by one. If an original
upper parameter is exactly 1, its shifted value N+1 cancels the new lower N+1,
so the rank need not increase. This includes the family with one exactly-unit
upper parameter relevant to the validation examples.

An exactly terminating upper parameter -m makes all n>m coefficients zero when
the lower parameters are admissible. That is a finite polynomial and needs no
infinite-tail algorithm. Exact singular lower parameters, including simultaneous
upper/lower nonpositive integers, must not be cancelled indiscriminately:
parameter-limit continuation and terminating-series conventions can differ.

## 5. Conservative pole budgets

For an affine parameter \(u=u_0+u_1\epsilon\), with
\(u_0\in\mathbb Z_{\le0}\) and \(u_1\ne0\), the valuation of \((u)_n\) is
one exactly when \(n\ge1-u_0\); otherwise it is zero. There is only one vanishing
factor in that Pochhammer product. Hence the generic coefficient valuation is

\[
v_n=\sum_i\mathbf1_{n\ge1-a_{i0},\,a_{i0}\le0}
    -\sum_j\mathbf1_{n\ge1-b_{j0},\,b_{j0}\le0}.
\]

An input-function pole bound is \(B=\max(0,-\min_n v_n)\); it can be obtained
by checking the finitely many threshold values and n=0. The simpler count of
nonpositive lower parameters is always a conservative upper bound in this
nonterminating affine setting. These bounds do not cover nonlinear regulator
parameters or exactly singular lower parameters.

For implementation, preserve a separate **intermediate shift budget**: a chain
with r divisions by epsilon needs its base through target+r, regardless of final
cancellations. External prefactors with pole bound B_ext require still deeper
orders. Positive shifts do not add pole budget. Do not subtract numerator zeros
from a local stage budget unless the cancellation occurs exactly before that
stage's truncation.

## Verification scope

`tests/test_integer_shift_math.py` checks coefficient identities, finite-prefix
reconstruction, positive-lower integral coefficients, and the missing-order
counterexample with exact SymPy arithmetic. These tests establish the algebra
for representative generic parameters. They do not certify a completed general
pFq translator, branch continuation or endpoint integration.

### Moving negative lower parameters in the tail construction

The same tail identity is valid for lower parameters \(-m+\beta\epsilon\)
with nonzero beta, for epsilon away from zero and away from other denominator
zeros. It therefore defines a meromorphic identity in a punctured neighborhood
of zero. Choose N greater than the magnitude of every nonpositive integer part
in both upper and lower lists; all tail parameters then have positive integer
parts. If \(C_N\) has valuation v, its regular tail must be expanded through
\(K+\max(0,-v)\) to recover the product through \(\epsilon^K\). When an
external factor has a pole bound B, use \(K+B+\max(0,-v)\), or a sharper
valuation of the exact combined prefactor if rigorously known. All finite-prefix
terms must also be expanded with the external prefactor, not discarded.

For example, upper parameters \((-1+\epsilon,1+2\epsilon,2+3\epsilon)\)
and lower parameters \((-2+4\epsilon,-3+5\epsilon)\), with N=4, give a
\(C_N\) with one simple pole: one zero numerator factor versus two zero
denominator factors. The shifted rank-4/rank-3 tail is regular. The tests verify
the exact coefficient identity and explicitly show a one-order-short expansion
corrupting the requested highest coefficient. An identically-zero C_N from an
exact terminating upper parameter must be handled as a finite polynomial first;
do not infer a finite valuation for the zero expression.

## Generic contiguous engine

`contiguous.py` implements algebraic parameter shifts directly from the
hypergeometric differential operator, rather than increasing the rank of a
negative-parameter tail. It plans the extra regulator depth before applying any
operator. Rational Laurent coefficient convolution retains apparent poles
until their cancellation is checked. An identically vanishing shift divisor
uses one formal slope perturbation; its removable limit must be established
exactly, otherwise the calculation fails. No numerical limit is substituted.

`canonical.py` selects ordinary nested sums, colored sums, or rationalized
half-integer differential systems. This is a Python-native implementation of
mathematical reduction ideas described by Huber and Maître,
[hep-ph/0507094](https://arxiv.org/abs/hep-ph/0507094) and
[0708.2443](https://arxiv.org/abs/0708.2443), not copied Mathematica source.

## Symbolic regulator slopes

Affine parameters may use exact polynomial slopes in independent symbols,
for example `1+a*epsilon` or `(a*b+1)*epsilon`. Slopes must not depend on the
hypergeometric argument. For a moving nonpositive integer lower parameter
`-m+b*epsilon`, an explicit nonzero assumption on `b` is required.

At every summation index, its Pochhammer product has at most one zero-base
factor, namely `b*epsilon`. Consequently the product D of these lower slopes
(including multiplicities) clears every genuine parameter denominator in the
regulator coefficients. Hypra checks that each HPL-basis coefficient becomes
polynomial in all slope symbols after multiplication by D. Other contiguous
operator denominators must cancel, so equality/zero specializations of regular
parameters do not inherit spurious singularities. This is a mathematical guard,
not an assumption that arbitrary slope differences are nonzero.
