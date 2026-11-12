"""Exact real HPL transform t=(1-z)/(1+z), with endpoint-derived constants."""

from functools import lru_cache

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import integrate_from_zero, shuffle_words
from .hpl_alternating import hpl_alternating_at_one
from .hpl_transform import _validate
from .symbolic import HPL


@lru_cache(maxsize=512)
def _mobius(word, z):
    if not word:
        return s.S.One
    one = HPL((-1,), z) - HPL((0,), z) - s.log(2)
    if all(letter == 1 for letter in word):
        return one ** len(word) / s.factorial(len(word))
    if word[0] == 1:
        counts = shuffle_words((1,), word[1:])
        diagonal = counts.pop(word)
        return s.expand(
            (
                one * _mobius(word[1:], z)
                - sum(count * _mobius(other, z) for other, count in counts.items())
            )
            / diagonal
        )
    # These two leading letters have convergent endpoints at one. Their
    # pulled-back kernels are regular at zero; logarithmic tails are integrable.
    kernel = -1 / (1 + z)
    if word[0] == 0:
        kernel -= 1 / (1 - z)
    return hpl_alternating_at_one(word) + integrate_from_zero(kernel * _mobius(word[1:], z), z)


def mobius_hpl(word, z, *, domain=None):
    """Return H_word((1-z)/(1+z)) in HPLs at z on the real interval 0<z<1.

    Every expanded {-1,0,1} word through weight six is supported, subject to
    the exact integrator's resource budget. Endpoint constants are preserved
    as exact alternating multiple-zeta values when not reduced. The transform
    maps the real interval to itself; no complex continuation is claimed.
    """
    word = _validate(word, z, domain)
    if len(word) > 6:
        raise UnsupportedExpressionError("Mobius transformation supports weight at most six")
    return _mobius(word, z)
