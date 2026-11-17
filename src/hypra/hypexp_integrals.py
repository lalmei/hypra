"""Exact logarithmic moments for the native HypExp interface."""

from functools import cache, lru_cache
from typing import Literal

import sympy as s

from ._typing import ExactInteger, ExactScalar
from .errors import UnsupportedExpressionError
from .hpl.algebra import integrate_from_zero, shuffle_expand


def _integers(*values: ExactInteger) -> tuple[int, ...]:
    if any(isinstance(value, bool) or not isinstance(value, (int, s.Integer)) for value in values):
        raise TypeError("Moment powers must be explicit integers")
    return tuple(map(int, values))


def beta_log_moment(n: ExactInteger, m: ExactInteger, p: ExactInteger) -> s.Expr:
    r"""Return the integral of u**p log(u)**n log(1-u)**m from zero to one.

    The origin requires p + m >= 0. Coefficients of a beta-function generating
    series give exact rational combinations of zeta values. Log weight through
    32 and polynomial degrees through 4096 are supported; these are algebra
    budgets, not restrictions of the beta-function identity.

    Args:
        n: Nonnegative power of the origin logarithm.
        m: Nonnegative power of the endpoint logarithm.
        p: Integer power of the integration variable.

    Returns:
        Exact value of the definite integral.

    Raises:
        TypeError: A power is not an explicit integer.
        ValueError: A logarithm power is negative.
        UnsupportedExpressionError: The integral diverges or exceeds the supported bounds.
    """
    n, m, p = _integers(n, m, p)
    if n < 0 or m < 0:
        raise ValueError("Logarithm powers must be nonnegative")
    if p + m < 0:
        raise UnsupportedExpressionError("The beta-log moment diverges at zero")
    if n + m > 32 or p > 4096:
        raise UnsupportedExpressionError("Beta-log moment exceeds the weight or power budget")
    return _beta_log_moment(n, m, p)


def _exp_bivariate(logarithm, n, m):
    """Exponentiate a rectangular exact Taylor jet with zero constant term.

    Euler's homogeneous-degree operator gives
    (i+j) F[i,j] = sum (r+t) L[r,t] F[i-r,j-t] for F = exp(L).
    This avoids differentiating large symbolic Gamma expressions.
    """
    coefficients = {(0, 0): s.S.One}
    for i in range(n + 1):
        for j in range(m + 1):
            if i == j == 0:
                continue
            value = sum(
                (
                    (r + t) * coefficient * coefficients[i - r, j - t]
                    for (r, t), coefficient in logarithm.items()
                    if r <= i and t <= j
                ),
                s.S.Zero,
            )
            coefficients[i, j] = s.expand(value / (i + j))
    return coefficients


@lru_cache(maxsize=256)
def _beta_log_moment(n, m, p):
    r"""Extract n! m! [a^n b^m] B(p+1+a,1+b), including removable poles.

    For p>=0 the logarithm of the beta function follows from the Taylor
    expansion of log Gamma at positive integers. For p=-1 write B(a,1+b)
    as E(a,b)/a, with E=Gamma(1+a)Gamma(1+b)/Gamma(1+a+b).

    For p=-r-1, r>=1, Gamma recurrence gives

        B(a-r,1+b) = -(1+b/a) F(a,b)/r,
        F = E/(1-a/r) prod_{j=1}^{r-1} (1+b/(a-j)).

    The residue in a is a polynomial in b of degree r. The convergence
    condition m>=r+1 kills it before coefficient extraction. Therefore these
    coefficients equal the convergent log moment, not a regularized divergent
    integral. The finite products enter through generalized harmonic numbers;
    no HPL tables, endpoint limits, or numerical fitting are used.
    """
    if m == 0:
        return s.S.NegativeOne**n * s.factorial(n) / s.Integer(p + 1) ** (n + 1)
    upper_n = n if p >= 0 else n + 1
    harmonic_base = p + 1 if p >= 0 else max(0, -p - 2)
    harmonics = {k: s.harmonic(harmonic_base, k) for k in range(1, upper_n + m + 1)}
    logarithm = {}
    for i in range(upper_n + 1):
        for j in range(m + 1):
            if i == j == 0:
                continue
            weight = i + j
            if j == 0:
                if p >= 0:
                    value = s.S.NegativeOne**i / (i * s.Integer(p + 1) ** i)
                elif p < -1:
                    value = 1 / (i * s.Integer(-p - 1) ** i)
                else:
                    value = s.S.Zero
            elif i == 0:
                sign = s.S.NegativeOne**j if p >= 0 else -1
                value = sign * harmonics[j] / j
            else:
                zeta_coefficient = s.S.NegativeOne ** (weight + 1) * s.zeta(weight)
                harmonic_sign = s.S.NegativeOne**weight if p >= 0 else -1
                value = (
                    s.binomial(weight, i)
                    * (zeta_coefficient + harmonic_sign * harmonics[weight])
                    / weight
                )
            if value != 0:
                logarithm[i, j] = value
    coefficients = _exp_bivariate(logarithm, upper_n, m)
    factorials = s.factorial(n) * s.factorial(m)
    if p >= 0:
        return s.expand(factorials * coefficients[n, m] / (p + 1))
    if p == -1:
        return s.expand(factorials * coefficients[n + 1, m])
    return s.expand(-factorials * (coefficients[n, m] + coefficients[n + 1, m - 1]) / (-p - 1))


def hypergeometric_log_moment(
    p: ExactInteger,
    n: ExactInteger,
    m: ExactInteger,
    k: ExactInteger,
    q: ExactInteger,
    z: ExactScalar,
    *,
    argument_domain: Literal["real", "principal-cut-plane"] = "real",
) -> s.Expr:
    r"""Integrate u**p log(u)**n log(1-u)**m log(1-z*u)**k / (z*u-1)**q.

    Integration runs from zero to one. Powers are nonnegative integers.
    Endpoints z=0 and z=1 are handled exactly when convergent. Symbolic z
    is interpreted on the real interval 0<z<1 unless it is provably negative.
    Select argument_domain="principal-cut-plane" for exact finite nonreal
    points, using the principal logarithm continued from z=0.

    Args:
        p: Power of the integration variable.
        n: Power of the origin logarithm.
        m: Power of the endpoint logarithm.
        k: Power of the moving-argument logarithm.
        q: Power of the moving denominator.
        z: Exact real kernel argument, or a symbol with the stated domain.
        argument_domain: Real kernels by default, or principal continuation
            to exact numeric points off the cut [1,infinity).

    Returns:
        Exact value of the definite integral.

    Raises:
        TypeError: A power is not an explicit integer.
        ValueError: A power is negative.
        UnsupportedExpressionError: The argument, powers, or convergence are unsupported.
    """
    p, n, m, k, q = _integers(p, n, m, k, q)
    if argument_domain not in ("real", "principal-cut-plane"):
        raise ValueError("argument_domain must be real or principal-cut-plane")
    if min(p, n, m, k, q) < 0:
        raise ValueError("Integral powers must be nonnegative")
    if p > 64 or q > 32 or n + m + k + int(q > 0) > 5:
        raise UnsupportedExpressionError("Logarithmic integral exceeds the weight or power budget")

    z = s.sympify(z, strict=True)
    if not isinstance(z, s.Expr) or z.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("The integral argument must be exact and finite")
    if argument_domain == "principal-cut-plane" and (z.free_symbols or z.is_finite is not True):
        raise UnsupportedExpressionError(
            "Principal-cut-plane moments require an exact finite numeric argument"
        )
    if z.is_zero is True:
        return s.S.Zero if k else s.S.NegativeOne**q * beta_log_moment(n, m, p)
    if z == 1:
        if n < q:
            raise UnsupportedExpressionError("The logarithmic integral diverges at one")
        # Reflect u about 1/2, then expand the polynomial (1-u)**p.
        return s.Add(
            *(
                s.S.NegativeOne ** (q + j) * s.binomial(p, j) * beta_log_moment(m + k, n, j - q)
                for j in range(p + 1)
            )
        )
    if k == q == 0:
        return beta_log_moment(n, m, p)
    if z.is_negative is True:
        # With v=1-u and y=z/(z-1), 1-z*u=(1-z)*(1-y*v).
        # Both factors are positive, so the logarithms add on the real branch.
        y = s.cancel(z / (z - 1))
        logarithm = s.log(1 - z)
        return (
            s.Add(
                *(
                    (-1) ** degree
                    * s.binomial(p, degree)
                    * s.binomial(k, weight)
                    * logarithm ** (k - weight)
                    * hypergeometric_log_moment(degree, m, n, weight, q, y)
                    for degree in range(p + 1)
                    for weight in range(k + 1)
                )
            )
            / (1 - z) ** q
        )
    if not (argument_domain == "principal-cut-plane" and z.is_real is False):
        if z.is_real is False or z.is_positive is False or (1 - z).is_positive is False:
            raise UnsupportedExpressionError("The moving-kernel formula requires 0 < z < 1")
        if z.is_number and (z.is_positive is not True or (1 - z).is_positive is not True):
            raise UnsupportedExpressionError("Could not establish 0 < z < 1")
    # For nonreal z, the radial segment t=z*u crosses neither the principal
    # log cut nor t=1. The real-interval germs below therefore continue to
    # these points by analyticity. log(t)-log(z)=log(u) along this segment;
    # log(1-z*u) is the principal branch since its segment avoids (-infinity,0].
    if m:
        coordinate = s.Dummy("argument", positive=True)
        moment = _mixed_log_moment(p, n, m, k, q, coordinate)
        return s.S.NegativeOne**q * moment.xreplace({coordinate: z})
    # Set t=z*u so all t-dependent poles are at 0 and 1.
    t = s.Dummy("t", positive=True)
    integrand = t**p * (s.log(t) - s.log(z)) ** n * s.log(1 - t) ** k / (t - 1) ** q
    return integrate_from_zero(integrand, t).subs(t, z) / z ** (p + 1)


def _mixed_log_moment(p, n, m, k, q, z):
    """Build the exact origin germ with positive denominator (1-z*u)**q."""

    def theta(expression):
        return z * s.diff(expression, z)

    @cache
    def ordinary(p, n, m, k):
        if k == 0:
            return beta_log_moment(n, m, p)
        if p:
            # Integrate the derivative of
            # u**(p+a)*(1-u)**(1+b)*(1-z*u)**c and differentiate
            # in (a,b,c). The endpoints vanish near (0,0,0).
            previous = ordinary(p - 1, n, m, k)
            result = p * previous + (z - 1) * s.diff(previous, z)
            if n:
                result += n * (ordinary(p - 1, n - 1, m, k) - ordinary(p, n - 1, m, k))
            if m:
                result -= m * ordinary(p, n, m - 1, k)
            result -= k * ordinary(p, n, m, k - 1)
            return shuffle_expand(result / (p + 1), z)
        # Differentiate the Euler equation for
        # B(1+a,1+b) * 2F1(-c,1+a;2+a+b;z).
        # Its zeroth-parameter operator is (1-z)*theta*(theta+1).
        previous = ordinary(0, n, m, k - 1)
        rhs = -k * z * (theta(previous) + previous)
        if n:
            rhs += n * (z - 1) * theta(ordinary(0, n - 1, m, k))
            rhs -= n * k * z * ordinary(0, n - 1, m, k - 1)
        if m:
            rhs -= m * theta(ordinary(0, n, m - 1, k))
        first = integrate_from_zero(rhs / (z * (1 - z)), z)
        return integrate_from_zero(first, z) / z

    @cache
    def moment(k, q):
        if q == 0:
            return ordinary(p, n, m, k)
        if q == 1:
            # theta(log(1-z*u)**(k+1)) gives the difference of the
            # q=1 and q=0 integrands, multiplied by -(k+1).
            result = moment(k, 0) - theta(moment(k + 1, 0)) / (k + 1)
        else:
            # Differentiation raises the denominator power. The log
            # derivative contributes only at one lower logarithm power.
            previous = moment(k, q - 1)
            result = previous + theta(previous) / (q - 1)
            if k:
                result += k * (moment(k - 1, q) - moment(k - 1, q - 1)) / (q - 1)
        return shuffle_expand(result, z)

    return moment(k, q)
