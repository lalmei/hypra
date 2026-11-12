"""Canonical half-integer hypergeometric epsilon systems.

Independent differential-equation implementation, not upstream package code.
The solution is the analytic germ at z=0, continued on the principal sheet.
The rationalizing variables use principal square roots. In case A, the path
for 0<z<1 is the positive imaginary t axis; the integrals are analytic germs,
not an assertion that a real-path integral may be substituted arbitrarily.
"""
from __future__ import annotations

import sympy as s

from .hpl_algebra import integrate_from_zero


def canonical_half_series(upper_base, upper_slopes, lower_base, lower_slopes,
                          z, epsilon, order):
    """Return a canonical expansion, or None for an unsupported base pattern.

    Slopes may be exact symbolic expressions independent of z and epsilon.
    Only pF(p-1), p>=2 is covered. The returned polynomial includes order.
    Integration/order failures propagate, never returning a partial series.
    """
    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("order must be a nonnegative integer")
    a, b = tuple(map(s.sympify, upper_base)), tuple(map(s.sympify, lower_base))
    alpha, beta = tuple(map(s.sympify, upper_slopes)), tuple(map(s.sympify, lower_slopes))
    p = len(a)
    if len(alpha) != p or len(beta) != len(b):
        raise ValueError("Each parameter needs one slope")
    if p < 2 or len(b) != p - 1:
        return None
    if any(x.has(epsilon, s.Float, s.oo, -s.oo, s.zoo, s.nan) or (isinstance(z, s.Symbol) and x.has(z))
           for x in alpha + beta):
        raise ValueError("Slopes must be exact and independent of argument and regulator")
    half = s.Rational(1, 2)
    if all(x == half for x in a+b):
        return _all_half_series(alpha, beta, z, epsilon, order)
    if all(x == 0 for x in a) and b.count(half) == 1 and all(x in (half, 1) for x in b):
        kind = 'A'
    elif a.count(half) == 1 and all(x in (half, 0) for x in a) and all(x == 1 for x in b):
        kind = 'B'
    elif p == 2 and sorted(a) == [0, half] and b == (half,):
        kind = 'C'
    else:
        return None
    t = s.Dummy('t', positive=True)
    x = s.Dummy('theta')
    if kind == 'A':
        zt = t**2 / (t**2 - 1)
        theta_factor = t * (1-t**2) / 2
        argument = s.I * s.sqrt(z/(1-z))
    elif kind == 'B':
        zt = 4*t/(1+t)**2
        theta_factor = t*(1+t)/(1-t)
        argument = (1-s.sqrt(1-z))/(1+s.sqrt(1-z))
    else:
        zt = t**2
        theta_factor = t/2
        argument = s.sqrt(z)
    operator = s.Poly(s.expand(x*s.prod(x+v-1+w*epsilon for v,w in zip(b,beta))
                              - zt*s.prod(x+v+w*epsilon for v,w in zip(a,alpha))), epsilon)

    def theta(value):
        return s.expand(theta_factor * s.diff(value, t))

    def apply(poly, value):
        if poly == 0 or value == 0:
            return s.S.Zero
        coefficients = s.Poly(poly, x)
        powers = [value]
        for _ in range(coefficients.degree()):
            powers.append(theta(powers[-1]))
        return s.expand(sum(c*powers[int(m[0])] for m,c in coefficients.terms()))

    def invert(value):
        value = s.expand(value)
        if value == 0:
            return s.S.Zero
        if kind == 'A':
            result = t * integrate_from_zero(2*value/t**2, t)
        elif kind == 'B':
            result = (1+t)/(1-t) * integrate_from_zero(value/t, t)
        else:
            result = t/(1-t**2) * integrate_from_zero(2*value/t**2, t)
        for _ in range(p-1):
            result = integrate_from_zero(s.expand(result/theta_factor), t)
        return s.expand(result)

    coefficients = [s.S.One]
    for n in range(1, order+1):
        forcing = -sum(apply(operator.nth(k), coefficients[n-k])
                       for k in range(1, min(n, operator.degree())+1))
        coefficients.append(invert(forcing))
    return sum(value.subs(t, argument)*epsilon**n for n,value in enumerate(coefficients))


def _all_half_series(alpha, beta, z, epsilon, order):
    """Gauge-normalized all-half pF(p-1) germ, using only ordinary HPL letters.

    L0=(theta-1/2)^(p-1) B in this order, B=(1-z)theta-z/2.
    Gauge F=(1-z)^(-1/2)Y, z=t²/(t²-1): conjugated theta
    is t(1-t²)/2 d/dt -t²/2. Conjugated B=t/2 d/dt, and
    conjugated (theta-1/2)=t(1-t²)/2 d/dt -(1+t²)/2.
    Both inverse kernels have only 0,±1 poles. Origin analyticity excludes
    the fractional-power homogeneous solutions of the latter operator.
    """
    t, x = s.Dummy('t', positive=True), s.Dummy('theta')
    half = s.Rational(1, 2)
    zt = t**2/(t**2-1)
    operator = s.Poly(s.expand(x*s.prod(x-half+b*epsilon for b in beta)
                              -zt*s.prod(x+half+a*epsilon for a in alpha)), epsilon)

    def apply(poly, value):
        if poly == 0 or value == 0:
            return s.S.Zero
        polynomial = s.Poly(poly, x)
        powers = [value]
        for _ in range(polynomial.degree()):
            last = powers[-1]
            powers.append(s.expand(t*(1-t**2)*s.diff(last,t)/2-t**2*last/2))
        return s.expand(sum(c*powers[int(m[0])] for m,c in polynomial.terms()))

    coefficients = [s.S.One]
    for n in range(1, order+1):
        forcing = -sum(apply(operator.nth(k), coefficients[n-k])
                       for k in range(1,min(n,operator.degree())+1))
        for _ in range(len(beta)):
            forcing = t/(1-t**2)*integrate_from_zero(s.expand(2*forcing/t**2),t)
        coefficients.append(integrate_from_zero(s.expand(2*forcing/t),t))
    argument = s.I*s.sqrt(z/(1-z))
    return sum(c.subs(t,argument)*epsilon**n for n,c in enumerate(coefficients))/s.sqrt(1-z)
