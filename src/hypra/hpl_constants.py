"""Fast convergent Euler-sum constants from exact real Mobius identities.

No direct, slowly convergent endpoint summation and no implicit branch choice.
Each unreduced endpoint is evaluated using HPLs at 1/3 and 1/2, recursively
requiring constants of strictly smaller weight. Precision checks are evidence,
not interval certification of all errors.
"""

from dataclasses import dataclass
from functools import lru_cache

import mpmath as mp
import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_alternating import hpl_alternating_at_one, reduce_euler_sum
from .hpl_endpoints import reduce_mzv
from .hpl_mobius import mobius_hpl
from .hpl_numeric import HPLNumericalConvergenceError, hpl_evalf
from .symbolic import HPL, MZV, EulerSum


@dataclass(frozen=True)
class ConstantNumericResult:
    """Approximate real constant; stability_gap is not an interval error bound."""

    value: object
    stability_gap: object
    dps: int
    working_dps: int
    method: str = 'real-Mobius-weight-recursion'
    convention: str = 'strict outer-index cutoff; no regularization'


def _word(indices, colors):
    word, cumulative, factor = [], 1, 1
    for index, color in zip(indices, colors, strict=True):
        cumulative *= int(color)
        word.extend([0] * (int(index) - 1))
        word.append(cumulative)
        factor *= cumulative
    return tuple(word), factor


@lru_cache(maxsize=256)
def _remainder(word):
    z = s.Symbol('constant_z', positive=True)
    remainder = s.expand(mobius_hpl(word, z, domain='0<z<1')
                         - hpl_alternating_at_one(word))
    for atom in remainder.atoms(EulerSum, MZV):
        if sum(atom.indices) >= len(word):
            raise UnsupportedExpressionError('Endpoint recursion did not strictly lower weight')
    return remainder.subs(z, s.Rational(1, 2))


def _calculate(constant, digits, max_terms):
    ctx = mp.mp.clone()
    ctx.dps = digits

    def as_float(value):
        return s.Float(ctx.nstr(value, n=digits), digits)

    @lru_cache(maxsize=512)
    def evaluate(atom):
        if atom == 1:
            return s.S.One
        weight = sum(atom.indices)
        if weight > 6:
            raise UnsupportedExpressionError('Numerical endpoint constants support weight at most six')
        if isinstance(atom, MZV):
            reduced = reduce_mzv(atom)
            if reduced == atom:
                return evaluate(EulerSum(atom.indices, (1,) * len(atom.indices)))
        else:
            reduced = reduce_euler_sum(atom)
        # Unreduced same-weight MZVs convert back to this EulerSum; evaluate
        # those by the endpoint recursion rather than creating a cycle.
        if reduced != atom and not reduced.has(EulerSum, MZV):
            return reduced.evalf(digits)
        colors = tuple(atom.colors) if isinstance(atom, EulerSum) else (1,) * len(atom.indices)
        word, factor = _word(atom.indices, colors)
        expression = _remainder(word)
        constants = expression.atoms(EulerSum, MZV)
        if any(sum(child.indices) >= weight for child in constants):
            raise UnsupportedExpressionError('Endpoint recursion did not strictly lower weight')
        expression = expression.xreplace({child: evaluate(child) for child in constants})
        expression = expression.xreplace({
            h: as_float(hpl_evalf(tuple(h.word), h.argument, dps=digits, max_terms=max_terms))
            for h in expression.atoms(HPL)
        })
        residual = expression.evalf(digits)
        if residual.free_symbols or residual.is_real is not True or residual.is_finite is not True:
            raise UnsupportedExpressionError('Endpoint recursion left an unevaluated numeric expression')
        value = hpl_evalf(word, s.Rational(1, 3), dps=digits, max_terms=max_terms)
        return as_float(factor * (ctx.convert(value) - ctx.mpf(str(residual))))

    result = evaluate(constant).evalf(digits)
    return ctx.mpf(str(result))


def _run(constant, dps, max_terms):
    if type(dps) is not int or not 2 <= dps <= 900:
        raise ValueError('dps must be an integer between 2 and 900')
    if type(max_terms) is not int or not 1 <= max_terms <= 1000000:
        raise ValueError('max_terms must be an integer between 1 and 1000000')
    if constant != 1 and sum(constant.indices) > 6:
        raise UnsupportedExpressionError('Numerical endpoint constants support weight at most six')
    first = _calculate(constant, dps + 20, max_terms)
    second = _calculate(constant, dps + 40, max_terms)
    ctx = mp.mp.clone()
    ctx.dps = dps + 40
    gap = abs(ctx.convert(first) - ctx.convert(second))
    if gap > ctx.power(10, -dps - 4) * abs(second):
        raise HPLNumericalConvergenceError('Endpoint constant unstable at increased precision')
    return ConstantNumericResult(second, gap, dps, dps + 40)


def evaluate_euler_sum(indices, colors, *, dps=50, max_terms=100000):
    """Numerically evaluate a convergent strict colored sum of weight <=6.

    The two-list input uses EulerSum's explicit sigma**n convention. Each
    recursive step must lower weight. HPL tail targets and repeated guarded
    precision checks are not a certified bound on the final constant's error.
    """
    return _run(EulerSum(indices, colors), dps, max_terms)


def evaluate_mzv(indices, *, dps=50, max_terms=100000):
    """Numerically evaluate a convergent strict multiple zeta of weight <=6."""
    return _run(MZV(indices), dps, max_terms)
