"""Exact binary HPL square and fractional argument maps."""

from functools import lru_cache

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_reduce import _extract
from .hpl_transform import _validate
from .symbolic import HPL


def _check(word, z, domain):
    word = _validate(word, z, domain)
    if any(letter not in (0, 1) for letter in word):
        raise UnsupportedExpressionError("This argument map supports only binary HPL words")
    if len(word) > 10:
        raise UnsupportedExpressionError("Argument map supports weight at most ten")
    return word


@lru_cache(maxsize=1024)
def _pullback(word, mode):
    """Expand a tensor word under the exact pullback of its one-forms."""
    if not word:
        return (((), s.S.One),)
    images = (
        {0: ((0, 2),), 1: ((1, 1), (-1, -1))}
        if mode == "square"
        else {0: ((0, 1), (1, 1)), 1: ((1, -1),)}
    )
    return tuple(
        ((letter, *tail), factor * coefficient)
        for letter, factor in images[word[0]]
        for tail, coefficient in _pullback(word[1:], mode)
    )


def _expression(word, z, mode):
    return s.Add(
        *(coefficient * (HPL(w, z) if w else s.S.One) for w, coefficient in _pullback(word, mode))
    )


def square_hpl(word, z, *, domain=None):
    """Express H_word(z**2) at z for binary words, real 0<z<1.

    The alphabet expands to {-1,0,1}; all trailing-zero constants follow from
    log(z**2)=2*log(z) on this declared interval.
    """
    return _expression(_check(word, z, domain), z, "square")


def fractional_hpl(word, z, *, boundary, domain=None):
    """Express H_word(z/(z-1) +/- i0), with explicit upper/lower boundary.

    Original letters must be binary. Mixed trailing zeros are extracted by
    shuffle before imposing log(z/(z-1))=H_0(z)+H_1(z)+/-i*pi.
    Origin-analytic words have identical boundary values on this negative ray.
    """
    word = _check(word, z, domain)
    if boundary not in ("upper", "lower"):
        raise UnsupportedExpressionError("Fractional boundary must be 'upper' or 'lower'")
    if not word:
        return s.S.One
    if word[-1] != 0:
        return _expression(word, z, "fractional")
    sign = 1 if boundary == "upper" else -1
    logarithm = HPL((0,), z) + HPL((1,), z) + sign * s.I * s.pi
    if all(letter == 0 for letter in word):
        return logarithm ** len(word) / s.factorial(len(word))
    extracted = _extract(word, z, False, True)
    return s.expand(
        extracted.xreplace(
            {
                atom: fractional_hpl(atom.word, z, boundary=boundary, domain=domain)
                for atom in extracted.atoms(HPL)
            }
        )
    )
