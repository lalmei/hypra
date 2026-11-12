"""Exact Taylor coefficients at nonsingular real interior base points."""
import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import _word
from .symbolic import HPL


def hpl_taylor(word, center, variable, *, through=3):
    """Expand HPL(word, variable) about exact rational 0<center<1.

    Coefficients retain exact HPL constants at center. The kernel ODE yields a
    finite convolution for every coefficient; no numerical differentiation is
    used. The Taylor series converges at least for |variable-center| less than
    min(center,1-center). Includes the coefficient of degree ``through`` plus
    the explicit Order remainder. Weight <=24 and degree <=64 are supported.
    """
    word = _word(word)
    if not isinstance(variable, s.Symbol) or variable.is_commutative is not True:
        raise TypeError('variable must be a commutative SymPy Symbol')
    if type(center) is int:
        center = s.Integer(center)
    if not isinstance(center, s.Rational) or not 0 < center < 1:
        raise UnsupportedExpressionError('Taylor center must be an exact rational in (0,1)')
    if type(through) is not int or not 0 <= through <= 64:
        raise UnsupportedExpressionError('Taylor degree must be an integer from zero through 64')
    coefficients = [s.S.One] + [s.S.Zero]*through
    for offset in range(len(word)-1, -1, -1):
        letter = word[offset]
        kernel = [(1 if letter == 1 else (-1)**j) /
                  ((1-center if letter == 1 else center if letter == 0 else 1+center)**(j+1))
                  for j in range(through)]
        coefficients = [HPL(word[offset:], center)] + [
            s.expand(sum(kernel[j]*coefficients[n-1-j] for j in range(n))/n)
            for n in range(1, through+1)]
    delta = variable-center
    polynomial = sum(coefficient*delta**n for n,coefficient in enumerate(coefficients))
    return polynomial + s.Order(delta**(through+1), (variable, center))
