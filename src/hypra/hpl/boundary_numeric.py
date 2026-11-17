"""Numerical real-cut boundary values from exact symbolic continuation maps."""

from dataclasses import dataclass

import mpmath as mp
import sympy as s

from .._typing import Boundary, ExactScalar, IntegerWord, RealRegion
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL, MZV, EulerSum
from .constants import evaluate_euler_sum, evaluate_mzv
from .numeric import HPLNumericalConvergenceError, MPComplex, MPReal, _word, evaluate_hpl
from .radial_continuation import _transport
from .real_continuation import analytic_continue_hpl


@dataclass(frozen=True)
class HPLBoundaryNumericResult:
    """Boundary value and non-certified precision diagnostics.

    component_tail_bounds are individual interior-HPL series truncation bounds,
    NOT a bound on the transformed result. For Taylor transport, tail_bound
    instead bounds propagated truncation. Neither includes roundoff error.
    """

    value: MPComplex
    stability_gap: MPReal
    component_tail_bounds: tuple[MPReal, ...]
    dps: int
    working_dps: int
    region: RealRegion
    boundary: Boundary
    method: str = "exact-real-boundary-map-to-interior"
    path: str = "principal boundary limit; no extra winding"
    tail_bound: MPReal | None = None


def _point(argument, region):
    if type(argument) is int:
        argument = s.Integer(argument)
    if (
        not isinstance(argument, s.Expr)
        or argument.free_symbols
        or argument.has(s.Float)
        or argument.is_real is not True
        or argument.is_finite is not True
    ):
        raise UnsupportedExpressionError(
            "Boundary argument must be an exact finite real SymPy number"
        )
    conditions = {
        "0<x<1": (argument.is_positive, (1 - argument).is_positive),
        "-1<x<0": ((-argument).is_positive, (argument + 1).is_positive),
        "x>1": ((argument - 1).is_positive,),
        "x<-1": ((-argument - 1).is_positive,),
    }
    if region not in conditions or any(test is not True for test in conditions[region]):
        raise UnsupportedExpressionError(
            "Exact argument does not satisfy the declared open real region"
        )
    return argument


def _evaluate(expression, digits, max_terms):
    ctx = mp.mp.clone()
    ctx.dps = digits
    replacements, tails = {}, []
    for atom in expression.atoms(HPL):
        result = evaluate_hpl(tuple(atom.word), atom.argument, dps=digits, max_terms=max_terms)
        replacements[atom] = s.Float(ctx.nstr(result.value, n=digits), digits)
        tails.append(result.tail_bound)
    for atom in expression.atoms(EulerSum, MZV):
        result = (
            evaluate_euler_sum(atom.indices, atom.colors, dps=digits, max_terms=max_terms)
            if isinstance(atom, EulerSum)
            else evaluate_mzv(atom.indices, dps=digits, max_terms=max_terms)
        )
        replacements[atom] = s.Float(ctx.nstr(result.value, n=digits), digits)
    result = expression.xreplace(replacements).evalf(digits)
    if (
        result.has(HPL, EulerSum, MZV)
        or result.free_symbols
        or result.is_number is not True
        or result.is_finite is not True
    ):
        raise UnsupportedExpressionError(
            "Boundary transformation left an unevaluated numerical expression"
        )
    real, imag = result.as_real_imag()
    return ctx.mpc(str(real.evalf(digits)), str(imag.evalf(digits))), tuple(tails)


def _boundary_transport(word, argument, region, boundary, dps, max_terms):
    """Continue along a finite detour ending exactly on the specified cut.

    Starting at 1/4, the three-segment path lies in the chosen open half
    plane until it reaches the target. It avoids 0, +/-1 and is homotopic
    to the requested principal boundary approach. The target has no small
    imaginary displacement. Interior real targets use a straight segment.
    """
    initial = s.Rational(1, 4)
    if region == "0<x<1":
        targets = [argument]
    else:
        height = s.I / 2 if boundary == "upper" else -s.I / 2
        targets = [initial + height, argument + height, argument]
    first = _transport(word, initial, targets, dps, dps + 25, 4096, 4096, max_terms)
    second = _transport(word, initial, targets, dps, dps + 45, 4096, 4096, max_terms)
    ctx = mp.mp.clone()
    ctx.dps = dps + 45
    gap = abs(ctx.mpc(first[0]) - ctx.mpc(second[0]))
    if gap > ctx.power(10, -dps - 4) * max(ctx.one, abs(second[0])):
        raise HPLNumericalConvergenceError(
            "Boundary transport unstable under increased working precision"
        )
    return HPLBoundaryNumericResult(
        ctx.mpc(second[0]),
        gap,
        (),
        dps,
        dps + 45,
        region,
        boundary,
        method="triangular-suffix-ODE-Taylor-boundary",
        path=(
            "real interior straight path"
            if region == "0<x<1"
            else "principal boundary limit via finite half-plane detour; no extra winding"
        ),
        tail_bound=second[1],
    )


def evaluate_hpl_boundary(
    word: IntegerWord,
    argument: ExactScalar,
    *,
    region: RealRegion,
    boundary: Boundary,
    dps: int = 50,
    max_terms: int = 100000,
) -> HPLBoundaryNumericResult:
    """Evaluate H_word(x +/- i0) by exact maps, never a tiny imaginary offset.

    Both region and boundary are mandatory. Regions: '0<x<1', '-1<x<0',
    'x>1', 'x<-1'. Boundaries: 'upper' (+i0), 'lower' (-i0). The real
    argument must satisfy the declared open region exactly. Singular endpoints
    and extra winding paths are excluded. All ternary words through weight
    24 are supported. Above weight six, a finite
    half-plane detour transports the defining ODE to the exact real target.
    Taylor transport also supplies a fallback if interior series do not
    converge within their term budget.

    The reported component bounds cover interior series tails individually,
    not propagation/cancellation/roundoff in the final boundary value. For
    Taylor transport, tail_bound gives the propagated truncation bound.
    Two guarded precision calculations must agree; this is not interval proof.

    Args:
        word: Expanded HPL letters.
        argument: Exact finite real point in the declared region.
        region: Open real interval containing the argument.
        boundary: Side of approach to the real axis.
        dps: Requested decimal precision, from 2 to 800.
        max_terms: Maximum terms per interior HPL series.

    Returns:
        Boundary value with stability and component truncation diagnostics.
    """
    if type(dps) is not int or not 2 <= dps <= 800:
        raise ValueError("dps must be an integer between 2 and 800")
    if type(max_terms) is not int or not 1 <= max_terms <= 1000000:
        raise ValueError("max_terms must be an integer between 1 and 1000000")
    argument = _point(argument, region)
    word = _word(word)
    if boundary not in ("upper", "lower"):
        raise UnsupportedExpressionError("Continuation boundary must be 'upper' or 'lower'")
    if len(word) > 6 and region in ("x>1", "x<-1"):
        return _boundary_transport(word, argument, region, boundary, dps, max_terms)
    x = s.Symbol("boundary_argument", real=True)
    expression = analytic_continue_hpl(word, x, region=region, boundary=boundary).subs(x, argument)
    try:
        first, _ = _evaluate(expression, dps + 20, max_terms)
        second, tails = _evaluate(expression, dps + 40, max_terms)
    except HPLNumericalConvergenceError:
        return _boundary_transport(word, argument, region, boundary, dps, max_terms)
    ctx = mp.mp.clone()
    ctx.dps = dps + 40
    gap = abs(ctx.convert(first) - ctx.convert(second))
    if gap > ctx.power(10, -dps - 4) * max(ctx.one, abs(second)):
        raise HPLNumericalConvergenceError(
            "Boundary value unstable under increased working precision"
        )
    return HPLBoundaryNumericResult(second, gap, tails, dps, dps + 40, region, boundary)
