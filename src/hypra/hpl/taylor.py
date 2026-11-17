"""Exact Taylor coefficients at nonsingular explicitly branched base points."""

from typing import Literal

import sympy as s

from .._typing import ExactScalar, IntegerWord
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL
from .algebra import _word


def hpl_taylor(
    word: IntegerWord,
    center: ExactScalar,
    variable: s.Symbol,
    *,
    through: int = 3,
    domain: Literal["real", "principal-cut-plane"] = "real",
) -> s.Expr:
    """Expand HPL(word, variable) about an exact nonsingular center.

    Coefficients retain exact HPL constants at center. The kernel ODE yields a
    finite convolution for every coefficient; no numerical differentiation is
    used. The default center is real with 0<center<1; exact algebraic and
    transcendental numbers are permitted. Explicit principal-cut-plane mode
    also permits exact nonreal centers. The Taylor series continues the local
    germ at center and converges inside the disk excluding 0, +/-1; it does
    not reset the branch when that disk crosses an artificial logarithm cut.
    Includes degree ``through`` and an Order remainder. Weight <=24 and
    degree <=64 are supported.

    Args:
        word: Expanded ternary word of weight at most 24.
        center: Exact finite numerical expansion point in the chosen domain.
        variable: Symbolic argument of the Taylor series.
        through: Inclusive Taylor degree, from zero to 64.
        domain: Real interior by default, or explicit principal cut plane.
    """
    word = _word(word)
    if not isinstance(variable, s.Symbol) or variable.is_commutative is not True:
        raise TypeError("variable must be a commutative SymPy Symbol")
    if type(center) is int:
        center = s.Integer(center)
    if domain not in ("real", "principal-cut-plane"):
        raise UnsupportedExpressionError("Taylor domain must be real or principal-cut-plane")
    if (
        not isinstance(center, s.Expr)
        or center.free_symbols
        or center.has(s.Float)
        or center.is_finite is not True
    ):
        raise UnsupportedExpressionError("Taylor center must be an exact finite number")
    interior = (
        center.is_real is True and center.is_positive is True and (1 - center).is_positive is True
    )
    if not interior and not (domain == "principal-cut-plane" and center.is_real is False):
        raise UnsupportedExpressionError(
            "Taylor center must be interior real or explicitly off-cut"
        )
    if type(through) is not int or not 0 <= through <= 64:
        raise UnsupportedExpressionError("Taylor degree must be an integer from zero through 64")
    coefficients = [s.S.One] + [s.S.Zero] * through
    for offset in range(len(word) - 1, -1, -1):
        letter = word[offset]
        kernel = [
            (1 if letter == 1 else (-1) ** j)
            / ((1 - center if letter == 1 else center if letter == 0 else 1 + center) ** (j + 1))
            for j in range(through)
        ]
        coefficients = [HPL(word[offset:], center)] + [
            s.expand(sum(kernel[j] * coefficients[n - 1 - j] for j in range(n)) / n)
            for n in range(1, through + 1)
        ]
    delta = variable - center
    polynomial = sum(coefficient * delta**n for n, coefficient in enumerate(coefficients))
    return polynomial + s.Order(delta ** (through + 1), (variable, center))
