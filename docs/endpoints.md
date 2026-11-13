# Endpoint and terminating-series preprocessing

`rewrite_endpoints(expr, epsilon, assumptions)` runs before Laurent-order
planning. It uses the convergence conditions and identities below.

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
hypergeometric endpoint. Higher-rank unit endpoints use the reductions below.

## Sources and implementation provenance

Gauss's identity: NIST DLMF 15.4.20,
https://dlmf.nist.gov/15.4.E20 . The defining generalized hypergeometric
series is DLMF 16.2.1, https://dlmf.nist.gov/16.2.E1 .

The broader expansion context is Huber and Maitre, *HypExp, a Mathematica
package for expanding hypergeometric functions around integer-valued
parameters*, hep-ph/0507094, and *HypExp 2*, arXiv:0708.2443. These endpoint
routines are independent implementations of the mathematical identities;
no upstream package source or lookup tables are copied.

## Higher-rank integer and half-integer endpoints

For supported affine integer and half-integer parameter families with
strictly positive base excess, the native engine expands at a dummy
interior argument and evaluates the limit of each complete Laurent
coefficient. Integer kernels use binary HPL reflection. Half-integer
kernels use the local coordinate `t = sqrt(1-z)` and power/log expansions.
The matched and single-half-upper coordinates tend to one. The
single-half-lower and all-half systems use an imaginary coordinate that
tends to infinity. The upper-side inversion identity is continued from
large positive arguments through the upper quadrant to obtain its
power/log expansion. These methods apply within the resource budgets of the interior kernels.
The limit is extracted after combining all
potentially divergent terms. The near-one engine supports HPL weight
through eight, and the inversion engine supports weight through six.

### Zero base excess on a convergent regulator side

Let $c_n$ be the defining coefficient at $z=1$ and let
$\Delta=\sum b_j-\sum a_i$. The sum of
$(n+1)c_{n+1}-n c_n$ vanishes when $\Re\Delta>0$, since
$n c_n=O(n^{-\Re\Delta})$. Adding $\Delta c_n$ to this
difference gives $c_n R(n)/\prod_j(n+b_j)$, where

\[
R(n)=\prod_i(n+a_i)+(\Delta-n)\prod_j(n+b_j).
\]

The two highest powers cancel. Expanding $R(n)$ in falling factorials
and shifting the summation index expresses the original series as a
linear combination of same-rank series with excess $\Delta+1$, divided
by $\Delta$. This rewrite runs before Laurent planning so the new pole
and any external poles are budgeted together. The zero-base-excess
rewrite requires assumptions proving $\Re\Delta>0$ on a punctured
regulator side. Negative base excess and identically zero excess
remain unsupported.

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
