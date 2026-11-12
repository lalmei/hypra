"""Generic exact contiguous shifts of pF(p-1) canonical regulator expansions.

Let L=A(theta)-z B(theta), A=theta*prod(theta+b-1), B=prod(theta+a).
The inverse shifts follow by polynomial division, not an inverse differential
operator with an unspecified integration constant. At a coefficient z**n they
multiply by (a-1)/(n+a-1) or b/(n+b), respectively.

All epsilon orders are propagated coefficientwise. An identically zero shift
divisor is regularized by ONE formal slope perturbation with distinct nonzero
integer weights, cancelled exactly and taken to zero algebraically. No numeric
parameter limit or fixed singular lower parameter is accepted.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import sympy as s

from .errors import IncompleteExpansionError, UnsupportedExpressionError
from .symbolic import HPL

_MAX_SHIFTS = 24
_MAX_PARAMETERS = 8
_MAX_DEPTH = 24
_MAX_TERMS = 50000


@dataclass(frozen=True)
class ShiftStep:
    kind: str
    index: int
    upper_bases: tuple
    lower_bases: tuple
    order_loss: int


@dataclass(frozen=True)
class ContiguousPlan:
    steps: tuple[ShiftStep, ...]
    extra_order: int
    perturbation_used: bool
    pole_bound: int


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise UnsupportedExpressionError("Contiguous bases must be exact rationals")
    return s.Rational(value)


def _parameters(upper_pairs, lower_pairs, canonical_upper_bases, canonical_lower_bases):
    from .nested_sums import _slope

    try:
        upper = tuple((_rational(base), _slope(slope)) for base, slope in upper_pairs)
        lower = tuple((_rational(base), _slope(slope)) for base, slope in lower_pairs)
        canonical_upper = tuple(map(_rational, canonical_upper_bases))
        canonical_lower = tuple(map(_rational, canonical_lower_bases))
    except (TypeError, ValueError) as exc:
        if isinstance(exc, UnsupportedExpressionError):
            raise
        raise UnsupportedExpressionError("Parameters must be (base, slope) pairs") from exc
    if not 2 <= len(upper) <= _MAX_PARAMETERS or len(upper) != len(lower) + 1:
        raise UnsupportedExpressionError("Contiguous engine requires 2 <= p <= 8 and pF(p-1)")
    if len(canonical_upper) != len(upper) or len(canonical_lower) != len(lower):
        raise UnsupportedExpressionError("Canonical and target parameter counts must match")
    if any(
        not (target - base).is_Integer
        for (target, _), base in zip(upper + lower, canonical_upper + canonical_lower, strict=True)
    ):
        raise UnsupportedExpressionError("Contiguous shifts must be integral")
    if any(base.is_Integer and base <= 0 and slope == 0 for base, slope in lower):
        raise UnsupportedExpressionError("Fixed singular lower parameter")
    if any(base.is_Integer and base <= 0 for base in canonical_lower):
        raise UnsupportedExpressionError("Canonical lower bases must be regular at epsilon=0")
    return upper, lower, canonical_upper, canonical_lower


def _divisor_factors(kind, index, upper_bases, lower_bases, upper_slopes, lower_slopes):
    if kind == "upper_up":
        return ((upper_bases[index], upper_slopes[index]),)
    if kind == "upper_down":
        return tuple(
            (base - upper_bases[index], slope - upper_slopes[index])
            for base, slope in zip(lower_bases, lower_slopes, strict=True)
        )
    if kind == "lower_down":
        return ((lower_bases[index] - 1, lower_slopes[index]),)
    return tuple(
        (base - lower_bases[index], slope - lower_slopes[index])
        for base, slope in zip(upper_bases, upper_slopes, strict=True)
    )


def plan_contiguous(
    upper_pairs, lower_pairs, canonical_upper_bases, canonical_lower_bases
) -> ContiguousPlan:
    """Plan shifts and their exact canonical-order budget without building a series."""
    upper, lower, upper_bases, lower_bases = _parameters(
        upper_pairs, lower_pairs, canonical_upper_bases, canonical_lower_bases
    )
    upper_bases, lower_bases = list(upper_bases), list(lower_bases)
    upper_slopes = tuple(slope for _, slope in upper)
    lower_slopes = tuple(slope for _, slope in lower)
    total_shifts = sum(
        abs(int(target - base))
        for (target, _), base in zip(
            upper + lower, tuple(upper_bases) + tuple(lower_bases), strict=True
        )
    )
    if total_shifts > _MAX_SHIFTS:
        raise UnsupportedExpressionError(f"Contiguous shift count exceeds {_MAX_SHIFTS}")
    steps = []
    perturbation = False
    for family, targets, bases in (("upper", upper, upper_bases), ("lower", lower, lower_bases)):
        for index, (target, _) in enumerate(targets):
            while bases[index] != target:
                direction = 1 if bases[index] < target else -1
                kind = family + ("_up" if direction == 1 else "_down")
                factors = _divisor_factors(
                    kind, index, upper_bases, lower_bases, upper_slopes, lower_slopes
                )
                loss = sum(1 for base, _ in factors if base == 0)
                perturbation |= any(base == 0 and slope == 0 for base, slope in factors)
                steps.append(ShiftStep(kind, index, tuple(upper_bases), tuple(lower_bases), loss))
                bases[index] += direction
    extra_order = sum(step.order_loss for step in steps)
    if extra_order > _MAX_DEPTH:
        raise UnsupportedExpressionError(f"Contiguous extra-order budget exceeds {_MAX_DEPTH}")
    return ContiguousPlan(
        tuple(steps),
        extra_order,
        perturbation,
        sum(1 for base, _ in lower if base.is_Integer and base <= 0),
    )


def _operator(step, upper_slopes, lower_slopes, z, epsilon):
    x = s.Dummy("theta")
    upper = [
        base + slope * epsilon for base, slope in zip(step.upper_bases, upper_slopes, strict=True)
    ]
    lower = [
        base + slope * epsilon for base, slope in zip(step.lower_bases, lower_slopes, strict=True)
    ]
    if step.kind == "upper_up":
        a = upper[step.index]
        operator = (x + a) / a
    elif step.kind == "upper_down":
        a = upper[step.index]
        A = s.Poly(x * s.prod(x + b - 1 for b in lower), x, domain="EX")
        quotient, _ = s.div(A, s.Poly(x + a - 1, x, domain="EX"))
        other_upper = s.prod(x + b for index, b in enumerate(upper) if index != step.index)
        operator = (quotient.as_expr() - z * other_upper) / s.prod(b - a for b in lower)
    elif step.kind == "lower_down":
        b = lower[step.index]
        operator = (x + b - 1) / (b - 1)
    else:
        b = lower[step.index]
        B = s.Poly(s.prod(x + a for a in upper), x, domain="EX")
        quotient, _ = s.div(B, s.Poly(x + b, x, domain="EX"))
        other_lower = x * s.prod(x + a - 1 for index, a in enumerate(lower) if index != step.index)
        operator = b * (other_lower - z * quotient.as_expr()) / (z * s.prod(a - b for a in upper))
    return tuple(
        (degree[0], s.cancel(coefficient))
        for degree, coefficient in s.Poly(s.expand(operator), x, domain="EX").terms()
    )


def _rational_laurent(expression, epsilon, through):
    """Exact regulator Laurent division of an operator's rational coefficient."""
    expression = s.cancel(expression)
    if expression == 0:
        return {}
    if not expression.is_rational_function(epsilon):
        raise UnsupportedExpressionError(
            "Contiguous operator coefficient is not rational in epsilon"
        )
    numerator, denominator = s.fraction(expression)
    numerator, denominator = (
        s.Poly(numerator, epsilon, domain="EX"),
        s.Poly(denominator, epsilon, domain="EX"),
    )
    nv = min(power[0] for power, coefficient in numerator.terms() if coefficient != 0)
    dv = min(power[0] for power, coefficient in denominator.terms() if coefficient != 0)
    valuation = nv - dv
    coefficients = []
    for degree in range(max(0, through - valuation + 1)):
        previous = sum(
            denominator.nth(dv + j) * coefficients[degree - j] for j in range(1, degree + 1)
        )
        coefficients.append(s.cancel((numerator.nth(nv + degree) - previous) / denominator.nth(dv)))
    return {
        valuation + degree: coefficient
        for degree, coefficient in enumerate(coefficients)
        if coefficient != 0
    }


def _coefficients(expression, epsilon, through):
    if not isinstance(expression, s.Expr) or expression.has(
        s.Order, s.Float, s.nan, s.zoo, s.oo, -s.oo
    ):
        raise IncompleteExpansionError(
            "Canonical builder must return an exact finite Laurent polynomial"
        )
    result = {}
    for term in s.Add.make_args(s.expand(expression, log=False, power_base=False, power_exp=False)):
        coefficient, power = term.as_coeff_exponent(epsilon)
        if not power.is_Integer or coefficient.has(epsilon):
            raise IncompleteExpansionError("Canonical builder left unresolved regulator dependence")
        if power < 0:
            raise IncompleteExpansionError(
                "Canonical regular lower bases unexpectedly produced regulator poles"
            )
        if power <= through:
            result[int(power)] = result.get(int(power), s.S.Zero) + coefficient
    return {
        power: s.expand(coefficient) for power, coefficient in result.items() if coefficient != 0
    }


def _hpl_pieces(expression):
    """Collect HPL atoms before rational cancellation to avoid denominator swell."""
    expression = s.expand(expression, log=False, power_base=False, power_exp=False)
    atoms = sorted(expression.atoms(HPL), key=s.default_sort_key)
    return (
        s.collect(expression, atoms, evaluate=False, exact=True) if atoms else {s.S.One: expression}
    )


def _normalize_coefficient(expression):
    return s.Add(
        *(s.cancel(coefficient) * word for word, coefficient in _hpl_pieces(expression).items())
    )


def _finite_delta_limit(expression, delta):
    """Exact removable limit, with no numeric extrapolation or arbitrary path."""
    expression = s.cancel(expression)
    if expression == 0 or not expression.has(delta):
        return expression
    if expression.is_rational_function(delta) is not True:
        raise UnsupportedExpressionError("Slope limit is not an exact rational expression")
    numerator, denominator = s.fraction(expression)
    numerator, denominator = (
        s.Poly(numerator, delta, domain="EX"),
        s.Poly(denominator, delta, domain="EX"),
    )
    nterms = [(power[0], s.simplify(coefficient)) for power, coefficient in numerator.terms()]
    dterms = [(power[0], s.simplify(coefficient)) for power, coefficient in denominator.terms()]
    nterms = [(power, coefficient) for power, coefficient in nterms if coefficient != 0]
    dterms = [(power, coefficient) for power, coefficient in dterms if coefficient != 0]
    if not nterms:
        return s.S.Zero
    if not dterms:
        raise UnsupportedExpressionError(
            "Slope perturbation produced an identically zero denominator"
        )
    nv, nc = min(nterms)
    dv, dc = min(dterms)
    if nv < dv:
        raise UnsupportedExpressionError("Degenerate contiguous slope limit did not cancel exactly")
    return s.cancel(nc / dc) if nv == dv else s.S.Zero


def expand_contiguous(
    upper_pairs,
    lower_pairs,
    canonical_upper_bases,
    canonical_lower_bases,
    z: s.Symbol,
    epsilon: s.Symbol,
    through: int,
    canonical_builder: Callable,
    *,
    max_internal_order: int = _MAX_DEPTH,
    assumptions=s.true,
) -> s.Expr:
    """Shift a canonical pF(p-1) expansion to exact rational target parameters.

    Each parameter is ``(base, slope)`` for base+slope*epsilon. Canonical bases
    must differ by integers from the targets and have regular lower parameters.
    The callback is

      canonical_builder((upper_bases, lower_bases), (upper_slopes, lower_slopes),
                        z, epsilon, inclusive_depth)

    and must return its complete exact polynomial through that depth. Slopes
    are exact polynomials in independent formal symbols; a zero divisor may
    require one additional polynomial delta. Moving singular lower slopes must be
    provably nonzero under assumptions. Multiplication by their product must
    clear all final slope denominators; unexplained parameter poles fail.
    The result has no Order term; the caller owns external prefactors and final
    expansion provenance. ``plan_contiguous`` exposes the required extra depth.
    """
    if not isinstance(z, s.Symbol) or not isinstance(epsilon, s.Symbol) or z == epsilon:
        raise TypeError("Use distinct symbolic argument and regulator variables")
    if type(through) is not int or through < 0:
        raise ValueError("through must be a nonnegative inclusive integer order")
    if type(max_internal_order) is not int or not 0 <= max_internal_order <= _MAX_DEPTH:
        raise ValueError(f"max_internal_order must be between 0 and {_MAX_DEPTH}")
    upper, lower, canonical_upper, canonical_lower = _parameters(
        upper_pairs, lower_pairs, canonical_upper_bases, canonical_lower_bases
    )
    plan = plan_contiguous(upper, lower, canonical_upper, canonical_lower)
    depth = through + plan.extra_order
    if depth > max_internal_order:
        raise UnsupportedExpressionError(
            f"Canonical depth {depth} exceeds contiguous limit {max_internal_order}"
        )
    slope_symbols = set().union(*(slope.free_symbols for _, slope in upper + lower))
    if {z, epsilon} & slope_symbols:
        raise UnsupportedExpressionError("Slopes must be independent of argument and regulator")
    from .planning import _facts

    singular_slopes = [slope for base, slope in lower if base.is_Integer and base <= 0]
    for slope in singular_slopes:
        if s.ask(s.Q.nonzero(slope), _facts(assumptions)) is not True:
            raise UnsupportedExpressionError("Moving lower pole slope must be provably nonzero")
    source_denominator = s.prod(singular_slopes)
    delta_name = "_hypra_slope_delta"
    while delta_name in {symbol.name for symbol in slope_symbols | {z, epsilon}}:
        delta_name += "_"
    delta = s.Symbol(delta_name)
    slopes = tuple(slope for _, slope in upper + lower)
    if plan.perturbation_used:
        slopes = tuple(slope + (index + 1) * delta for index, slope in enumerate(slopes))
    upper_slopes, lower_slopes = slopes[: len(upper)], slopes[len(upper) :]
    expression = canonical_builder(
        (canonical_upper, canonical_lower), (upper_slopes, lower_slopes), z, epsilon, depth
    )
    state = _coefficients(expression, epsilon, depth)
    for step in plan.steps:
        next_depth = depth - step.order_loss
        if not state:
            return s.S.Zero
        lower_power = min(state)
        operator = _operator(step, upper_slopes, lower_slopes, z, epsilon)
        max_derivative = max(degree for degree, _ in operator)
        derivatives = {0: state}
        for degree in range(1, max_derivative + 1):
            derivatives[degree] = {
                power: z * s.diff(value, z) for power, value in derivatives[degree - 1].items()
            }
        next_state = {}
        for degree, rational in operator:
            rational_series = _rational_laurent(rational, epsilon, next_depth - lower_power)
            for left_power, left_coefficient in rational_series.items():
                for right_power, right_coefficient in derivatives[degree].items():
                    power = left_power + right_power
                    if power <= next_depth:
                        next_state[power] = (
                            next_state.get(power, s.S.Zero) + left_coefficient * right_coefficient
                        )
        state = {}
        for power, coefficient in next_state.items():
            coefficient = _normalize_coefficient(coefficient)
            if len(s.Add.make_args(coefficient)) > _MAX_TERMS:
                raise UnsupportedExpressionError(
                    "Contiguous coefficient exceeds the supported term budget"
                )
            if coefficient != 0:
                state[power] = coefficient
        depth = next_depth
    result = s.S.Zero
    for power, coefficient in state.items():
        coefficient = _normalize_coefficient(coefficient)
        if plan.perturbation_used:
            coefficient = s.Add(
                *(
                    _finite_delta_limit(value, delta) * word
                    for word, value in _hpl_pieces(coefficient).items()
                )
            )
        if slope_symbols:
            for value in _hpl_pieces(coefficient).values():
                if s.cancel(source_denominator * value).is_polynomial(*slope_symbols) is not True:
                    raise IncompleteExpansionError("Uncanceled symbolic slope denominator")
        if coefficient.has(delta, s.nan, s.zoo, s.oo, -s.oo, epsilon):
            raise IncompleteExpansionError(
                "Contiguous result retained invalid or unresolved parameters"
            )
        if power < -plan.pole_bound:
            if s.simplify(coefficient) != 0:
                raise IncompleteExpansionError(
                    "Apparent contiguous regulator poles failed to cancel"
                )
            continue
        if power <= through:
            result += coefficient * epsilon**power
    return result
