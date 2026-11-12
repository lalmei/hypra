# Endpoint and terminating-series preprocessing

`rewrite_endpoints(expr, epsilon, assumptions)` runs before Laurent-order
planning. It performs only the identities below. It neither evaluates a
divergent endpoint by continuation nor exchanges a limit and an integral.

## Finite terminating sums

For a fixed upper parameter -m, m a nonnegative integer, the defining series
terminates. With M the smallest such m,

\[
 {}_pF_q(\mathbf a;\mathbf b;z)
 =\sum_{n=0}^{M}\frac{\prod_i(a_i)_n}{\prod_j(b_j)_n}\frac{z^n}{n!}.
\]

This finite algebraic expression is valid at any finite argument, including
z=1, where its lower-parameter denominators are defined. A moving parameter
-m+alpha epsilon does **not** terminate. A moving lower parameter may give a
meromorphic regulator pole; retaining the rational finite sum exposes that
pole to the subsequent planner. Any fixed nonpositive integer lower
parameter is rejected, even if termination appears to occur before its
first zero: competing limiting/termination conventions are not chosen
implicitly. Symbolic termination bounds and degrees above 4096 are outside
the current implementation budget.

## Gauss's convergent sum

For a nonterminating 2F1 at z=1, Gauss's summation theorem states

\[
 {}_2F_1(a,b;c;1)
 =\frac{\Gamma(c)\Gamma(c-a-b)}{\Gamma(c-a)\Gamma(c-b)},
 \qquad \Re(c-a-b)>0.
\]

The lower parameter must not be a fixed nonpositive integer. Parameters are
required to be affine in the regulator. If the regulator-zero excess is a
positive rational, fixed finite slopes admit an open neighborhood where
its real part remains positive. If that base is zero, the implementation
requires assumptions proving positive real excess on the selected punctured
side. For example, 2F1(1,1;2-epsilon;1) is convergent for real epsilon<0 and
its exact Gamma expression simplifies to (epsilon-1)/epsilon. Its pole is
therefore included in the prefactor-order budget after rewriting.

A negative base excess, identically zero excess, or unproven convergence is
rejected. The Gamma quotient may have a meromorphic continuation beyond the
convergent domain, but that fact is not used to assign a value to a divergent
hypergeometric endpoint. Higher-rank unit endpoints are restricted to the strictly convergent
all-integer path described below. Endpoint rewrites make no claim about integrals of the
result or uniform endpoint estimates needed to exchange other operations.

## Sources and implementation provenance

Gauss's identity: NIST DLMF 15.4.20,
https://dlmf.nist.gov/15.4.E20 . The defining generalized hypergeometric
series is DLMF 16.2.1, https://dlmf.nist.gov/16.2.E1 .

The broader expansion context is Huber and Maitre, *HypExp, a Mathematica
package for expanding hypergeometric functions around integer-valued
parameters*, hep-ph/0507094, and *HypExp 2*, arXiv:0708.2443. These endpoint
routines are independent implementations of the mathematical identities;
no upstream package source or lookup tables are copied.

## Generic all-integer pF(p-1) unit endpoints

The native engine also accepts integer-base affine parameters with rational
slopes when the base excess

\[
 \Delta_0=\sum_j b_j(0)-\sum_i a_i(0)>0.
\]

This path first expands at a dummy real argument 0<z<1 using the existing
integer algorithm. Each complete Laurent coefficient is then transformed
with z=1-t and binary HPL reflection. The finite origin-limit routine combines
all rational and logarithmic terms before deciding whether a finite limit
exists. It never adds individually divergent endpoint values.

The implementation is bounded by the reflection engine: binary HPL words
of weight at most eight and supported rational coefficients. Unsupported word
alphabets, arguments, weights or unproved cancellations raise. This is not a
generic half-integer endpoint evaluator. Nonpositive base excess still fails,
except when another explicitly supported identity such as Gauss or exact
termination applies.

### Why the coefficient endpoint limit is justified

For pF(p-1), the large-n defining coefficient behaves as
`C(epsilon) n^(-1-Delta(epsilon))` times a uniformly bounded asymptotic factor.
Choose a sufficiently small complex regulator disk so
Re Delta(epsilon)>Delta_0/2. At moving nonpositive lower bases, multiply by a
conservative regulator power epsilon^B to remove the finitely many parameter
poles. The Gamma-ratio prefactor is then analytic and bounded on a smaller
closed disk; fixed singular lower parameters remain rejected. Uniform Gamma
ratio asymptotics bound the regularized n-th coefficient by
`C n^(-1-Delta_0/2)` for all sufficiently large n, and the finitely many earlier
coefficients are bounded separately. Thus the regularized series converges
normally in the regulator and uniformly as real z approaches one from below.
Cauchy's coefficient formula allows extraction of each fixed regulator
coefficient and passage to that endpoint. Dividing by epsilon^B recovers the
Laurent expansion. No new endpoint regulator pole arises from this strictly
positive excess domain; moving lower poles remain accounted for by planning.

Examples checked independently in `tests/test_hpl_limits.py`:

\[
 {}_3F_2(1,1,1+\epsilon;2,2+\epsilon;1)
 =\zeta(2)+\epsilon[\zeta(2)-\zeta(3)]
 +\epsilon^2[\zeta(4)-\zeta(3)]+O(\epsilon^3).
\]

Its defining coefficient is
`(1+epsilon)/[(n+1)(n+1+epsilon)]`, so the displayed sums are elementary zeta
series, independent of HPL reflection. A moving-lower-pole case gives

\[
 {}_3F_2(-1+\epsilon,1,1;\epsilon,4;1)
 =-\frac1{4\epsilon}+\frac{53}{48}+O(\epsilon).
\]

Here n=1 supplies the pole, while the n>=2 constant tail is the convergent
rational sum `sum 6/[(n-1)(n+1)(n+2)(n+3)]=7/48`. An external 1/epsilon
prefactor test checks that endpoint evaluation preserves deeper-order planning.
The unit-series convergence criterion is also summarized in DLMF16.2:
https://dlmf.nist.gov/16.2 .
