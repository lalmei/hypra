"""Exact origin-continuation HPL values at i, without fitted constant tables."""

from functools import lru_cache

import sympy as s

from .._typing import IntegerWord
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL, reduce_known_hpl
from .classical import classical_hpl
from .multiple_polylog import MultiplePolylog
from .notation import PMWord, pm_to_letters
from .numeric import _word
from .reduction import _extract


@lru_cache(maxsize=512)
def _segment(word, start, displacement):
    """Iterated integral of a word with no singular final letter at start."""
    indices, arguments = [], []
    zeros, previous = 0, displacement
    for letter in word:
        pole = letter - start
        if pole == 0:
            zeros += 1
        else:
            indices.append(zeros + 1)
            arguments.append(s.expand(s.cancel(previous / pole)))
            zeros, previous = 0, pole
    sign = (-1) ** (len(indices) - word.count(1))
    return sign * MultiplePolylog(s.Tuple(*indices), s.Tuple(*arguments))


@lru_cache(maxsize=512)
def _first_half(word):
    if not word:
        return s.S.One
    midpoint = s.I / 2
    z = s.Symbol("_hpl_midpoint")
    extracted = _extract(word, z, False, True)
    return extracted.xreplace(
        {
            atom: s.log(midpoint) if atom.word == (0,) else _segment(tuple(atom.word), 0, midpoint)
            for atom in extracted.atoms(HPL)
        }
    )


def _at_i_multiple_polylogs(word):
    # Compose the first half with the reversal of the segment i -> i/2.
    # Both series have cumulative argument products with modulus <= 1/2.
    return _first_half(word) + s.Add(
        *(
            (-1) ** cut * _segment(word[:cut][::-1], s.I, -s.I / 2) * _first_half(word[cut:])
            for cut in range(1, len(word) + 1)
        )
    )


@lru_cache(maxsize=1093)
def _ordinary_at_i(word: tuple[int, ...]) -> s.Expr:
    """Exact H_word(i) on continuation from the origin through the upper quadrant.

    All ternary words through weight six reduce to exact logarithms and
    classical or multiple polylogarithms. Pure repeated letters and simple
    classical words are also supported through weight 24.

    Args:
        word: Expanded ternary word evaluated at the positive imaginary unit.
    """
    if not word:
        return s.S.One
    z = s.Dummy("real_formula_argument", positive=True)
    if len(word) <= 2:
        # Obtain the proved real formula first, then continue it analytically.
        # Along z=i*t, 0<t<=1, every log/Li2 argument in these low-weight
        # formulas avoids its cut. This is not a complex-domain assertion
        # passed to classical_hpl's real-only public interface.
        formula = classical_hpl(HPL(word, z), z, domain="0<z<1")
    else:
        formula = reduce_known_hpl(HPL(word, z))
        if formula.has(HPL):
            if len(word) <= 6:
                return _at_i_multiple_polylogs(word)
            raise UnsupportedExpressionError(
                "At i, mixed words above weight six are not implemented"
            )
    value = formula.subs(z, s.I)
    # The even/odd terms of Li2(±i) give -zeta(2)/8 ± i*Catalan.
    # Li2(1/2) follows from the reflection identity at 1/2.
    return s.expand(
        value.xreplace(
            {
                s.polylog(2, s.I): -(s.pi**2) / 48 + s.I * s.Catalan,
                s.polylog(2, -s.I): -(s.pi**2) / 48 - s.I * s.Catalan,
                s.polylog(2, s.S.Half): s.pi**2 / 12 - s.log(2) ** 2 / 2,
            }
        )
    )


def hpl_at_i(word: IntegerWord) -> s.Expr:
    """Return the exact ordinary HPL at i through weight six.

    Repeated single letters and simple classical words extend to weight 24.
    """
    return _ordinary_at_i(_word(word))


def hpl_pm_at_i(word: PMWord) -> s.Expr:
    """Return the exact plus/minus HPL at i through weight six."""
    terms = pm_to_letters(word)
    if len(word) > 6:
        raise UnsupportedExpressionError("At i, plus/minus words support weight at most six")
    return s.Add(*(coefficient * _ordinary_at_i(letters) for letters, coefficient in terms.items()))
