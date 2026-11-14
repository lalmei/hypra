"""Dispatch a single hypergeometric atom to an exact native kernel."""

import sympy as s

from ..errors import UnsupportedExpressionError
from ..symbolic import HPL, reduce_known_hpl
from .planning import _positive, plan_expansion


def _inside_interval(z, assumptions):
    return _positive(z, assumptions) and _positive(1 - z, assumptions)


def _linear(expr, eps, base):
    """Return rational slope if expr=base+slope*eps, else None."""
    slope = s.expand(expr).coeff(eps)
    if slope.is_Rational and s.expand(expr - base - slope * eps) == 0:
        return slope
    return None


def _zero_base_twof1(b, c, z, eps, order):
    # f1=b H_1; fk= -c H_(0,word)+(b-c) H_(1,word).
    result = s.S.One
    words = {(1,): b}
    for k in range(1, order + 1):
        result += eps**k * sum(coef * HPL(word, z) for word, coef in words.items())
        following = {}
        for word, coef in words.items():
            for letter, factor in ((0, -c), (1, b - c)):
                value = coef * factor
                if value:
                    following[(letter, *word)] = value
        words = following
    return reduce_known_hpl(result)


def _li(k, z):
    return -s.log(1 - z) if k == 1 else s.polylog(k, z)


def _half_telescoping(b, z, eps, order):
    t = s.sqrt(z)

    def odd_sum(k):
        return s.Integer(2) ** (k - 1) * (_li(k, t) - _li(k, -t)) / t

    result = (s.log(1 + t) - s.log(1 - t)) / (2 * t)
    for k in range(1, order + 1):
        result += (-1) ** (k - 1) * b**k * (odd_sum(k) - odd_sum(k + 1) / 2) * eps**k
    return result


def _half_even(a, z, eps, order):
    t = s.sqrt(z)
    return sum(
        ((-2 * a) ** k / (2 * s.factorial(k))) * (s.log(1 - t) ** k + s.log(1 + t) ** k) * eps**k
        for k in range(order + 1)
    )


def expand_atom(h, eps, order, assumptions=s.true, argument_domain="real-interval"):
    """Exact coefficient recursion in the explicit real domain 0<z<1."""
    if h.argument == 0:
        return s.S.One
    z = h.argument
    if z == 1:
        from .endpoints import convergent_unit
        from .planning import laurent_coefficient
        from .unit_limit import unit_coefficient_limit

        if not convergent_unit(h, eps):
            raise UnsupportedExpressionError("Unit endpoint requires strictly positive base excess")
        if order == 0:
            from .endpoints import rewrite_endpoints

            base = h.subs(eps, 0)
            if not any(b.is_Integer and b <= 0 for b in base.bq) and any(
                a.is_Integer and a <= 0 for a in base.ap
            ):
                return rewrite_endpoints(base, eps, assumptions)
        regular_z = s.Dummy("unit_argument", positive=True)
        regular_h = s.hyper(h.ap, h.bq, regular_z, evaluate=False)
        regular = expand_atom(regular_h, eps, order, s.And(assumptions, regular_z < 1), "real-interval")
        poles = plan_expansion(h, eps, order, assumptions).pole_bound
        return s.Add(*(eps**k * unit_coefficient_limit(laurent_coefficient(regular, eps, k), regular_z)
                       for k in range(-poles, order+1)))
    if argument_domain == "principal-disk":
        # This increment admits exact numerical complex points only. Symbolic
        # disk constraints require a separate assumption proof, not a guess.
        if z.free_symbols or not _positive(1-s.Abs(z), s.true) or not (
            z.is_real is False or _positive(z, s.true)
        ):
            raise UnsupportedExpressionError(
                "Principal-disk expansion requires an exact numeric argument with "
                "0 < |z| < 1 off the negative real axis"
            )
    elif argument_domain == "principal-cut-plane":
        if z.free_symbols or z.is_finite is not True or not (
            z.is_real is False or _inside_interval(z, s.true)
        ):
            raise UnsupportedExpressionError(
                "Principal-cut-plane expansion requires an exact finite numeric "
                "point off (-infinity,0] and [1,infinity)"
            )
    elif not _inside_interval(z, assumptions):
        raise UnsupportedExpressionError(
            "Native hypergeometric expansion requires 0 < argument < 1"
        )
    if len(h.ap) == 1 and not h.bq:
        return s.series((1 - z) ** (-h.ap[0]), eps, 0, order + 1).removeO()
    # Preserve compact classical forms for two common half-integer subfamilies.
    if len(h.ap) == 2 and len(h.bq) == 1:
        upper, lower = tuple(h.ap), h.bq[0]
        if s.S.One in upper:
            other = upper[1] if upper[0] == 1 else upper[0]
            b = _linear(other, eps, s.Rational(1, 2))
            c = _linear(lower, eps, s.Rational(3, 2))
            if b is not None and b == c:
                return _half_telescoping(b, z, eps, order)
        if lower == s.Rational(1, 2):
            for first, second in (upper, upper[::-1]):
                a = _linear(first, eps, s.S.Zero)
                b = _linear(second, eps, s.Rational(1, 2))
                if a is not None and a == b:
                    return _half_even(a, z, eps, order)
    from .canonical import build_canonical, canonical_bases, use_contiguous
    from .contiguous import expand_contiguous
    from .matched_half import half_parameters

    parsed = half_parameters(h, eps)
    if parsed is not None and use_contiguous(*parsed):
        ap, bq = parsed
        ca, cb = canonical_bases(ap, bq)
        coordinate = s.Dummy("argument", positive=True)
        result = expand_contiguous(
            ap, bq, ca, cb, coordinate, eps, order, build_canonical, assumptions=assumptions
        )
        return reduce_known_hpl(result).xreplace({coordinate: z})
    from .integer import expand_integer_hyper

    integer_result = expand_integer_hyper(h, eps, order)
    if integer_result is not None:
        return integer_result
    from .matched_half import expand_matched_half

    matched_result = expand_matched_half(h, eps, order)
    if matched_result is not None:
        return matched_result
    if len(h.ap) != 2 or len(h.bq) != 1:
        from .half_system import canonical_half_series
        from .matched_half import half_parameters

        parsed = half_parameters(h, eps)
        if parsed is not None:
            ap, bq = parsed
            result = canonical_half_series(
                [a for a, _ in ap], [a for _, a in ap],
                [b for b, _ in bq], [b for _, b in bq], z, eps, order,
            )
            if result is not None:
                return result
        raise UnsupportedExpressionError(
            "Unsupported noninteger parameter pattern or shifts"
        )
    upper, lower = tuple(h.ap), h.bq[0]
    if s.S.One in upper:
        other = upper[1] if upper[0] == 1 else upper[0]
        b, c = _linear(other, eps, s.S.Zero), _linear(lower, eps, s.S.One)
        if b is not None and c is not None:
            return _zero_base_twof1(b, c, z, eps, order)
        b = _linear(other, eps, s.Rational(1, 2))
        c = _linear(lower, eps, s.Rational(3, 2))
        if b is not None and b == c:
            return _half_telescoping(b, z, eps, order)
    if lower == s.Rational(1, 2):
        for first, second in (upper, upper[::-1]):
            a = _linear(first, eps, s.S.Zero)
            b = _linear(second, eps, s.Rational(1, 2))
            if a is not None and a == b:
                return _half_even(a, z, eps, order)
    from .twof1_shift import expand_shifted_twof1

    shifted = expand_shifted_twof1(h, eps, order)
    if shifted is not None:
        return shifted
    from .half_system import canonical_half_series
    from .matched_half import half_parameters

    parsed = half_parameters(h, eps)
    if parsed is not None:
        ap, bq = parsed
        canonical = canonical_half_series(
            [a for a, _ in ap],
            [a for _, a in ap],
            [b for b, _ in bq],
            [b for _, b in bq],
            z,
            eps,
            order,
        )
        if canonical is not None:
            return reduce_known_hpl(canonical)
    raise UnsupportedExpressionError(
        "Native port does not yet support this parameter family; no external fallback was used"
    )
