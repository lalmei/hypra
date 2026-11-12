"""Numerical real-cut boundary values from exact symbolic continuation maps."""

from dataclasses import dataclass

import mpmath as mp
import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_constants import evaluate_euler_sum, evaluate_mzv
from .hpl_continuation import analytic_continue_hpl
from .hpl_numeric import HPLNumericalConvergenceError, _word, evaluate_hpl
from .symbolic import HPL, MZV, EulerSum


@dataclass(frozen=True)
class HPLBoundaryNumericResult:
    """Boundary value and non-certified precision diagnostics.

    component_tail_bounds are individual interior-HPL series truncation bounds,
    NOT a bound on the transformed result. Stability excludes no roundoff risk.
    """

    value: object
    stability_gap: object
    component_tail_bounds: tuple
    dps: int
    working_dps: int
    region: str
    boundary: str
    method: str = 'exact-real-boundary-map-to-interior'
    path: str = 'principal boundary limit; no extra winding'


def _point(argument, region):
    if type(argument) is int:
        argument = s.Integer(argument)
    if (not isinstance(argument, s.Expr) or argument.free_symbols
            or argument.has(s.Float) or argument.is_real is not True
            or argument.is_finite is not True):
        raise UnsupportedExpressionError('Boundary argument must be an exact finite real SymPy number')
    conditions = {
        '0<x<1': (argument.is_positive, (1-argument).is_positive),
        '-1<x<0': ((-argument).is_positive, (argument+1).is_positive),
        'x>1': ((argument-1).is_positive,),
        'x<-1': ((-argument-1).is_positive,),
    }
    if region not in conditions or any(test is not True for test in conditions[region]):
        raise UnsupportedExpressionError('Exact argument does not satisfy the declared open real region')
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
        result = (evaluate_euler_sum(atom.indices, atom.colors, dps=digits, max_terms=max_terms)
                  if isinstance(atom, EulerSum) else
                  evaluate_mzv(atom.indices, dps=digits, max_terms=max_terms))
        replacements[atom] = s.Float(ctx.nstr(result.value, n=digits), digits)
    result = expression.xreplace(replacements).evalf(digits)
    if (result.has(HPL, EulerSum, MZV) or result.free_symbols
            or result.is_number is not True or result.is_finite is not True):
        raise UnsupportedExpressionError('Boundary transformation left an unevaluated numerical expression')
    real, imag = result.as_real_imag()
    return ctx.mpc(str(real.evalf(digits)), str(imag.evalf(digits))), tuple(tails)


def evaluate_hpl_boundary(word, argument, *, region, boundary, dps=50, max_terms=100000):
    """Evaluate H_word(x +/- i0) by exact maps, never a tiny imaginary offset.

    Both region and boundary are mandatory. Regions: '0<x<1', '-1<x<0',
    'x>1', 'x<-1'. Boundaries: 'upper' (+i0), 'lower' (-i0). The real
    argument must satisfy the declared open region exactly. Singular endpoints
    and extra winding paths are excluded. Outside [-1,1], weight <=6;
    inside, the symbolic continuation's weight limit is 24.

    The reported component bounds cover interior series tails individually,
    not propagation/cancellation/roundoff in the final boundary value. Two
    guarded precision calculations must agree; this is not interval proof.
    """
    if type(dps) is not int or not 2 <= dps <= 800:
        raise ValueError('dps must be an integer between 2 and 800')
    if type(max_terms) is not int or not 1 <= max_terms <= 1000000:
        raise ValueError('max_terms must be an integer between 1 and 1000000')
    argument = _point(argument, region)
    word = _word(word)
    x = s.Symbol('boundary_argument', real=True)
    expression = analytic_continue_hpl(word, x, region=region, boundary=boundary).subs(x, argument)
    first, _ = _evaluate(expression, dps+20, max_terms)
    second, tails = _evaluate(expression, dps+40, max_terms)
    ctx = mp.mp.clone()
    ctx.dps = dps+40
    gap = abs(ctx.convert(first)-ctx.convert(second))
    if gap > ctx.power(10, -dps-4)*abs(second):
        raise HPLNumericalConvergenceError('Boundary value unstable under increased working precision')
    return HPLBoundaryNumericResult(second, gap, tails, dps, dps+40, region, boundary)
