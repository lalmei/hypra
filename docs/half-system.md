# Canonical half-integer differential systems

This is an independent implementation of the hypergeometric differential
equation, not a translation of upstream Mathematica source. It covers the
three patterns below for rank p >= 2, not every HypExp2 half-integer class.
All parameters are base + slope × epsilon. Slopes may be symbolic and must
be independent of the argument and epsilon.

Write theta = z d/dz and

\[
 L(\epsilon)=\theta\prod_{j=1}^{p-1}(\theta+b_j-1)
              -z\prod_{i=1}^{p}(\theta+a_i).
\]

For F = sum epsilon^n F_n and L = sum epsilon^k L_k, the analytic solution
has F_0=1, F_n(0)=0 for n>0, and

\[
 L_0F_n=-\sum_{k=1}^nL_kF_{n-k}.
\]

## Rationalized inverse operators

In all cases L_0 = B theta^(p-1). Apply B inverse first, then theta inverse
p-1 times. All displayed primitives have their lower-origin constant fixed
by `integrate_from_zero`; a divergent primitive is not accepted.

**A:** every upper base is zero; one lower base is 1/2, all others 1.

\[
 z=\frac{t^2}{t^2-1},\quad
 \theta=\frac{t(1-t^2)}2\partial_t,\quad
 B=(1-z)\theta-\frac12,
\]
\[
 B^{-1}g=t\int_0^t\frac{2g(u)}{u^2}\,du,\qquad
 \theta^{-1}f=\int_0^t\frac{2f(u)}{u(1-u^2)}\,du.
\]

**B:** one upper base is 1/2, all others zero; all lower bases are 1.
The rank restriction is essential: at rank one F_0 would not equal 1.

\[
 z=\frac{4t}{(1+t)^2},\quad
 \theta=\frac{t(1+t)}{1-t}\partial_t,\quad
 B=(1-z)\theta-\frac z2,
\]
\[
 B^{-1}g=\frac{1+t}{1-t}\int_0^t\frac{g(u)}u\,du,\qquad
 \theta^{-1}f=\int_0^t\frac{1-u}{u(1+u)}f(u)\,du.
\]

**C:** 2F1 upper bases (0,1/2), lower base 1/2.

\[
 z=t^2,\quad \theta=\frac t2\partial_t,\quad
 B=(1-z)\theta-\frac{1+z}2,
\]
\[
 B^{-1}g=\frac{t}{1-t^2}\int_0^t\frac{2g(u)}{u^2}\,du,\qquad
 \theta^{-1}f=\int_0^t\frac{2f(u)}u\,du.
\]

These formulas follow by substituting the displayed integrating factors in
Bf=g. All rational poles are in the HPL alphabet {0,1,-1}. The RHS in A/C
is even and O(t²), so integration at the origin is regular. In B it is O(t),
also regular. Successive theta inverses preserve zero constant terms.
The homogeneous B solutions in A/C introduce odd t terms, corresponding to
fractional powers of z; analyticity in z excludes them. In B a homogeneous
B term has nonzero origin limit and subsequent inversion introduces logs;
F_n(0)=0 and analytic regularity exclude it. The kernels of repeated theta
are constants/logarithms, likewise excluded. Thus these inverses select the
unique defining hypergeometric analytic germ.

## Continuation and branch scope

The integrator establishes primitive identities for 0<t<1 with its explicit
HPL convention. Because the complete expressions above are analytic germs
at zero after removable singularities cancel, these same identities hold
in a complex neighborhood by analytic continuation; there is no interchange
of an improper endpoint integral and regulator expansion here.

For 0<z<1, B uses t=(1-sqrt(1-z))/(1+sqrt(1-z)) and C uses t=sqrt(z), both
in (0,1). A uses t=i sqrt(z/(1-z)), on the positive imaginary axis. Continuing
the origin germ into the upper half-plane and along this axis crosses
neither HPL singularity +1 nor -1. The complete germ is regular at zero,
so its continuation needs no prescription for a divergent origin integral.
The principal HPL continuation along that path is therefore the intended
value. This does **not** assert that an arbitrary complex branch or a path
winding around a singularity is equivalent. z=1 limits and endpoint
integrations require a separate analysis and are not justified here.

## Verification

`tests/test_half_system.py` expands the returned HPL words independently as
Taylor series in t using their defining kernels. It compares exact rational
coefficients through epsilon³ and t⁴ against finite Pochhammer-series
coefficients under each rational z(t) map, including rank-three A/B cases,
zero and negative slopes. These checks are independent of the implemented
ODE inverse, but finite checks are not a proof of unrestricted rank/order
performance or full HypExp2 coverage. Symbolic slopes and missing operator
coefficients are tested separately.
