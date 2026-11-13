"""Exact logarithmic moments for the native HypExp interface."""

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import integrate_from_zero
from .hpl_limits import hpl_limit_at_one


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
    return hpl_limit_at_one(primitive, u)


def hypergeometric_log_moment(p, n, m, k, q, z):
    r"""Integrate u**p log(u)**n log(1-u)**m log(1-z*u)**k / (z*u-1)**q.

    Integration runs from zero to one. Powers are nonnegative integers.
    Endpoints z=0 and z=1 are handled exactly. On 0<z<1 the moving-kernel
    formula requires m=0, except when both k and q vanish. Symbolic z is
    interpreted on that interval.
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
    if z.is_real is False or z.is_positive is False or (1-z).is_positive is False:
        raise UnsupportedExpressionError("The moving-kernel formula requires 0 < z < 1")
    if z.is_number and (z.is_positive is not True or (1-z).is_positive is not True):
        raise UnsupportedExpressionError("Could not establish 0 < z < 1")
    if m:
        raise UnsupportedExpressionError("The interior formula requires m=0 for a moving kernel")
    # Set t=z*u so all t-dependent poles are at 0 and 1.
    t = s.Dummy("t", positive=True)
    integrand = t**p * (s.log(t)-s.log(z))**n * s.log(1-t)**k / (t-1)**q
    return integrate_from_zero(integrand, t).subs(t, z) / z**(p+1)
