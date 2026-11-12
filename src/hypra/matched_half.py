"""Matched half-integer parameter reductions in the rational variable sqrt(z)."""

from __future__ import annotations

import sympy as s

from .errors import UnsupportedExpressionError
from .integer import _truncate
from .nested_sums import _slope
from .planning import _exact_rational_valuation, laurent_coefficient
from .symbolic import reduce_known_hpl


def half_parameters(h, eps):
    rows = []
    for parameters in (h.ap, h.bq):
        row = []
        for parameter in parameters:
            try:
                p = s.Poly(parameter, eps)
            except s.PolynomialError:
                return None
            base, slope = p.nth(0), p.nth(1)
            if (
                p.degree() > 1
                or not base.is_Rational
                or base.q not in (1, 2)
            ):
                return None
            try:
                slope = _slope(slope, eps)
            except (TypeError, ValueError):
                return None
            row.append((base, slope))
        rows.append(row)
    return tuple(rows)


def expand_matched_half(h, eps, order):
    """Expand matched upper/lower half counts, with arbitrary integer offsets."""
    from .colored_sums import canonical_colored_hyper_series
    from .hpl_algebra import integrate_from_zero

    parsed = half_parameters(h, eps)
    if parsed is None:
        return None
    upper, lower = parsed
    nh = sum(base.q == 2 for base, _ in upper)
    if nh == 0 or nh != sum(base.q == 2 for base, _ in lower):
        return None
    if len(upper) != len(lower) + 1:
        raise UnsupportedExpressionError("Matched-half engine requires pF(p-1)")
    if any(base.is_Integer and base <= 0 and slope == 0 for base, slope in lower):
        raise UnsupportedExpressionError("Fixed singular lower parameter")
    poles = sum(1 for base, _ in lower if base.is_Integer and base <= 0)
    terminating = [
        -int(base) for base, slope in upper if base.is_Integer and base <= 0 and slope == 0
    ]
    if terminating:
        value = sum(
            h.argument**n
            / s.factorial(n)
            * s.prod(s.rf(a, n) for a in h.ap)
            / s.prod(s.rf(b, n) for b in h.bq)
            for n in range(min(terminating) + 1)
        )
        return _truncate(value, eps, order, poles)
    if any(base < 0 for base, _ in upper):
        shift = int(s.floor(-min(base for base, _ in upper + lower))) + 1
        coefficient, prefix = s.S.One, s.S.Zero
        for n in range(shift):
            if n:
                coefficient *= s.prod(a + n - 1 for a in h.ap) / (
                    n * s.prod(b + n - 1 for b in h.bq)
                )
            prefix += coefficient * h.argument**n
        coefficient *= s.prod(a + shift - 1 for a in h.ap) / (
            shift * s.prod(b + shift - 1 for b in h.bq)
        )
        valuation = _exact_rational_valuation(coefficient, eps, s.true)
        tail_h = s.hyper(
            [1, *(a + shift for a in h.ap)],
            [shift + 1, *(b + shift for b in h.bq)],
            h.argument,
            evaluate=False,
        )
        tail = expand_matched_half(tail_h, eps, order + max(0, -valuation))
        return _truncate(prefix + coefficient * h.argument**shift * tail, eps, order, poles)
    t = s.Dummy("t", positive=True)
    theta = lambda f: t * s.diff(f, t) / 2
    depth = order + poles
    base_upper = [
        (s.Rational(1, 2) if base.q == 2 else (0 if base == 0 else 1), slope)
        for base, slope in upper
    ]
    base_lower = [(s.Rational(1, 2) if base.q == 2 else 1, slope) for base, slope in lower]
    result = canonical_colored_hyper_series(base_upper, base_lower, t**2, eps, depth)
    for (base, slope), (initial, _) in zip(upper, base_upper, strict=True):
        for step in range(int(base - initial)):
            a = initial + step + slope * eps
            result = _truncate((theta(result) + a * result) / a, eps, depth, poles)
    for (base, slope), (initial, _) in zip(lower, base_lower, strict=True):
        for step in range(int(initial - base)):
            new = initial - step - 1
            c = new + slope * eps
            if new == 0:
                depth -= 1
            result = _truncate((theta(result) + c * result) / c, eps, depth, poles)
    for (base, slope), (initial, _) in zip(lower, base_lower, strict=True):
        for step in range(int(base - initial)):
            current = initial + step
            expfactor = sum(
                (2 * slope * eps * s.log(t)) ** k / s.factorial(k) for k in range(depth + poles + 1)
            )
            weighted = _truncate(t ** (2 * current - 1) * expfactor * result, eps, depth, poles)
            integrated = s.S.Zero
            for k in range(-poles, depth + 1):
                coeff = laurent_coefficient(weighted, eps, k)
                if coeff:
                    integrated += eps**k * integrate_from_zero(coeff, t)
            inverse = sum(
                (-2 * slope * eps * s.log(t)) ** k / s.factorial(k)
                for k in range(depth + poles + 1)
            )
            result = _truncate(
                2 * (current + slope * eps) * t ** (-2 * current) * inverse * integrated,
                eps,
                depth,
                poles,
            )
    return reduce_known_hpl(_truncate(result, eps, order, poles)).xreplace({t: s.sqrt(h.argument)})
