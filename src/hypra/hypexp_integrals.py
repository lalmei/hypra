"""Exact logarithmic moments for the native HypExp interface."""

from functools import cache

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl.algebra import _linear_terms, integrate_from_zero, shuffle_expand
from .hpl.regularization import hpl_endpoint_polynomial


def _integers(*values):
    if any(isinstance(value, bool) or not isinstance(value, (int, s.Integer))
           for value in values):
        raise TypeError("Moment powers must be explicit integers")
    return values


def beta_log_moment(n, m, p):
    r"""Return the integral of u**p log(u)**n log(1-u)**m from zero to one.

    The origin requires p + m >= 0. The result uses exact HPL endpoint
    constants. The implementation supports log weight through seven, negative
    powers through -32, and polynomial degrees through 64.
    """
    n, m, p = _integers(n, m, p)
    if n < 0 or m < 0:
        raise ValueError("Logarithm powers must be nonnegative")
    if p + m < 0:
        raise UnsupportedExpressionError("The beta-log moment diverges at zero")
    if n + m > 7 or not -32 <= p <= 64:
        raise UnsupportedExpressionError("Beta-log moment exceeds the weight or power budget")
    u = s.Dummy("u", positive=True)
    primitive = integrate_from_zero(u**p * s.log(u)**n * s.log(1-u)**m, u)
    # Integration by parts produces rational coefficients with poles only
    # at zero. At one, replacing each coefficient by its value leaves an
    # error O((1-u)*log(1-u)**weight), even for divergent HPL factors.
    T = s.Symbol("_beta_endpoint_log")
    endpoint = s.expand(sum(
        coefficient.subs(u, 1)*hpl_endpoint_polynomial(word, T)
        for word, coefficient in _linear_terms(primitive, u).items()
    ))
    if endpoint.has(T, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("Could not establish a finite beta-log moment")
    return endpoint


def hypergeometric_log_moment(p, n, m, k, q, z):
    r"""Integrate u**p log(u)**n log(1-u)**m log(1-z*u)**k / (z*u-1)**q.

    Integration runs from zero to one. Powers are nonnegative integers.
    Endpoints z=0 and z=1 are handled exactly when convergent. Symbolic z
    is interpreted on the real interval 0<z<1 unless it is provably negative.
    """
    p, n, m, k, q = _integers(p, n, m, k, q)
    if min(p, n, m, k, q) < 0:
        raise ValueError("Integral powers must be nonnegative")
    if p > 64 or q > 32 or n + m + k + int(q > 0) > 5:
        raise UnsupportedExpressionError("Logarithmic integral exceeds the weight or power budget")

    z = s.sympify(z, strict=True)
    if not isinstance(z, s.Expr) or z.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("The integral argument must be exact and finite")
    if z.is_zero is True:
        return s.S.Zero if k else s.S.NegativeOne**q * beta_log_moment(n, m, p)
    if z == 1:
        if n < q:
            raise UnsupportedExpressionError("The logarithmic integral diverges at one")
        # Reflect u about 1/2, then expand the polynomial (1-u)**p.
        return s.Add(*(s.S.NegativeOne**(q+j) * s.binomial(p, j)
                       * beta_log_moment(m+k, n, j-q) for j in range(p+1)))
    if k == q == 0:
        return beta_log_moment(n, m, p)
    if z.is_negative is True:
        # With v=1-u and y=z/(z-1), 1-z*u=(1-z)*(1-y*v).
        # Both factors are positive, so the logarithms add on the real branch.
        y = s.cancel(z/(z-1))
        logarithm = s.log(1-z)
        return s.Add(*(
            (-1)**degree * s.binomial(p, degree) * s.binomial(k, weight)
            * logarithm**(k-weight)
            * hypergeometric_log_moment(degree, m, n, weight, q, y)
            for degree in range(p+1)
            for weight in range(k+1)
        )) / (1-z)**q
    if z.is_real is False or z.is_positive is False or (1-z).is_positive is False:
        raise UnsupportedExpressionError("The moving-kernel formula requires 0 < z < 1")
    if z.is_number and (z.is_positive is not True or (1-z).is_positive is not True):
        raise UnsupportedExpressionError("Could not establish 0 < z < 1")
    if m:
        coordinate = s.Dummy("argument", positive=True)
        moment = _mixed_log_moment(p, n, m, k, q, coordinate)
        return s.S.NegativeOne**q * moment.xreplace({coordinate: z})
    # Set t=z*u so all t-dependent poles are at 0 and 1.
    t = s.Dummy("t", positive=True)
    integrand = t**p * (s.log(t)-s.log(z))**n * s.log(1-t)**k / (t-1)**q
    return integrate_from_zero(integrand, t).subs(t, z) / z**(p+1)


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
            previous = ordinary(p-1, n, m, k)
            result = p*previous + (z-1)*s.diff(previous, z)
            if n:
                result += n*(ordinary(p-1, n-1, m, k)-ordinary(p, n-1, m, k))
            if m:
                result -= m*ordinary(p, n, m-1, k)
            result -= k*ordinary(p, n, m, k-1)
            return shuffle_expand(result/(p+1), z)
        # Differentiate the Euler equation for
        # B(1+a,1+b) * 2F1(-c,1+a;2+a+b;z).
        # Its zeroth-parameter operator is (1-z)*theta*(theta+1).
        previous = ordinary(0, n, m, k-1)
        rhs = -k*z*(theta(previous) + previous)
        if n:
            rhs += n*(z-1)*theta(ordinary(0, n-1, m, k))
            rhs -= n*k*z*ordinary(0, n-1, m, k-1)
        if m:
            rhs -= m*theta(ordinary(0, n, m-1, k))
        first = integrate_from_zero(rhs/(z*(1-z)), z)
        return integrate_from_zero(first, z)/z

    @cache
    def moment(k, q):
        if q == 0:
            return ordinary(p, n, m, k)
        if q == 1:
            # theta(log(1-z*u)**(k+1)) gives the difference of the
            # q=1 and q=0 integrands, multiplied by -(k+1).
            result = moment(k, 0) - theta(moment(k+1, 0))/(k+1)
        else:
            # Differentiation raises the denominator power. The log
            # derivative contributes only at one lower logarithm power.
            previous = moment(k, q-1)
            result = previous + theta(previous)/(q-1)
            if k:
                result += k*(moment(k-1, q)-moment(k-1, q-1))/(q-1)
        return shuffle_expand(result, z)

    return moment(k, q)
