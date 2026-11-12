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
hypergeometric endpoint. Nonterminating higher-rank unit endpoints currently
remain unsupported. Endpoint rewrites make no claim about integrals of the
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
