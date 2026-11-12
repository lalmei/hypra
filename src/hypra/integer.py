"""Native integer-parameter reduction using exact contiguous operators.

Independent implementation of standard coefficient identities, not copied
HypExp source. Canonical nested sums are expanded into HPLs before shifts.
"""

from __future__ import annotations

import sympy as s

from .errors import UnsupportedExpressionError
from .planning import laurent_coefficient
from .symbolic import reduce_known_hpl


def integer_parameters(h, eps):
    """Return exact integer bases/rational slopes, or None for a noninteger family."""
    result = []
    for parameters in (h.ap, h.bq):
        row = []
        for parameter in parameters:
            try:
                poly = s.Poly(parameter, eps)
            except s.PolynomialError:
                return None
            if poly.degree() > 1 or not poly.nth(0).is_Integer or not poly.nth(1).is_Rational:
                return None
            row.append((int(poly.nth(0)), poly.nth(1)))
        result.append(row)
    return tuple(result)


def _truncate(expr, eps, order, pole_bound=0):
    return s.expand(
        s.series(expr * eps**pole_bound, eps, 0, order + pole_bound + 1).removeO() / eps**pole_bound
    )


def _terminating(h, eps, order, upper, lower):
    degrees = [-base for base, slope in upper if base <= 0 and slope == 0]
    if not degrees:
        return None
    nmax = min(degrees)
    result = s.S.Zero
    for n in range(nmax + 1):
        term = h.argument**n / s.factorial(n)
        for p in h.ap:
            term *= s.rf(p, n)
        for p in h.bq:
            term /= s.rf(p, n)
        result += term
    poles = sum(1 for base, _ in lower if base <= 0)
    return _truncate(result, eps, order, poles)


def expansion_depth(h, eps, order):
    parsed = integer_parameters(h, eps)
    if parsed is None:
        return order
    _, lower = parsed
    return order + sum(1 for base, _ in lower if base <= 0)


def expand_integer_hyper(h, eps, order):
    """Expand integer bases with nonnegative upper offsets and any lower offset.

    Positive lower shifts use convergent integrals at zero. Lower nonpositive
    offsets require nonzero regulator slope and enough extra canonical order.
    Fixed negative upper parameters are handled as exact terminating series.
    Moving negative upper parameters use exact finite-prefix/tail reindexing,
    never a divergent inverse-upper integral.
    """
    from .nested_sums import canonical_hyper_series

    ap, bq = list(h.ap), list(h.bq)
    for b in tuple(bq):
        if b in ap and not (b.is_Integer and b <= 0):
            ap.remove(b)
            bq.remove(b)
    if len(ap) != len(h.ap):
        h = s.hyper(ap, bq, h.argument, evaluate=False)
    parsed = integer_parameters(h, eps)
    if parsed is None:
        return None
    upper, lower = parsed
    if any(base <= 0 and slope == 0 for base, slope in lower):
        raise UnsupportedExpressionError("Fixed singular lower parameter")
    terminating = _terminating(h, eps, order, upper, lower)
    if terminating is not None:
        return terminating
    if any(base < 0 for base, _ in upper):
        # Exact finite-prefix/tail reindexing avoids divergent inverse-upper
        # integrals. Shift all bases positive, allowing the canonical engine.
        from .planning import _exact_rational_valuation

        shift = 1 - min(base for base, _ in upper + lower)
        prefix = s.S.Zero
        coefficient = s.S.One
        for n in range(shift):
            if n:
                coefficient *= s.prod(a + n - 1 for a in h.ap)
                coefficient /= n * s.prod(b + n - 1 for b in h.bq)
            prefix += coefficient * h.argument**n
        coefficient *= s.prod(a + shift - 1 for a in h.ap)
        coefficient /= shift * s.prod(b + shift - 1 for b in h.bq)
        valuation = _exact_rational_valuation(coefficient, eps, s.true)
        shifted = s.hyper(
            [1, *(a + shift for a in h.ap)],
            [shift + 1, *(b + shift for b in h.bq)],
            h.argument,
            evaluate=False,
        )
        tail = expand_integer_hyper(shifted, eps, order + max(0, -valuation))
        return _truncate(
            prefix + coefficient * h.argument**shift * tail,
            eps,
            order,
            sum(1 for base, _ in lower if base <= 0),
        )
    if len(upper) != len(lower) + 1:
        raise UnsupportedExpressionError("Integer nested-sum engine requires pF(p-1)")
    if any(base <= 0 and slope == 0 for base, slope in lower):
        raise UnsupportedExpressionError("Fixed singular lower parameter")
    z = s.Dummy("z", positive=True)
    depth = expansion_depth(h, eps, order)
    canonical_upper = [(0 if base == 0 else 1, slope) for base, slope in upper]
    result = canonical_hyper_series(canonical_upper, [slope for _, slope in lower], z, eps, depth)

    # Raising an upper parameter a -> a+1: (theta+a)/a.
    for base, slope in upper:
        for current in range(1, base):
            a = current + slope * eps
            result = _truncate((z * s.diff(result, z) + a * result) / a, eps, depth)

    # Lowering c -> c-1: (theta+c-1)/(c-1). Each zero crossing
    # spends one additional base-series order; no O term is silently discarded.
    for base, slope in lower:
        for new in range(0, base - 1, -1):
            c = new + slope * eps
            if new == 0:
                depth -= 1
            result = _truncate(
                (z * s.diff(result, z) + c * result) / c,
                eps,
                depth,
                sum(1 for b, _ in lower if b <= 0),
            )

    # Raising positive c: c*z^(-c) int_0^z t^(c-1) F(t) dt.
    # Every coefficient is integrable because current >=1 and F(0)=1.
    for base, slope in lower:
        for current in range(1, base):
            from .hpl_algebra import integrate_from_zero

            c = current + slope * eps
            pole_bound = sum(1 for b, _ in lower if b <= 0)
            exponential = sum(
                (slope * eps * s.log(z)) ** k / s.factorial(k)
                for k in range(depth + pole_bound + 1)
            )
            weighted = _truncate(
                z ** (current - 1) * exponential * result,
                eps,
                depth,
                sum(1 for b, _ in lower if b <= 0),
            )
            integrated = s.S.Zero
            lower_order = -sum(1 for b, _ in lower if b <= 0)
            for k in range(lower_order, depth + 1):
                coeff = laurent_coefficient(weighted, eps, k)
                if coeff:
                    integrated += eps**k * integrate_from_zero(coeff, z)
            inverse = sum(
                (-slope * eps * s.log(z)) ** k / s.factorial(k)
                for k in range(depth - lower_order + 1)
            )
            result = _truncate(c * z ** (-current) * inverse * integrated, eps, depth, -lower_order)
    result = _truncate(result, eps, order, sum(1 for b, _ in lower if b <= 0))
    return reduce_known_hpl(result).xreplace({z: h.argument})
