"""Exact shuffle extraction of endpoint logarithms, without branch changes."""
from __future__ import annotations

from functools import lru_cache

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import _word, shuffle_words
from .symbolic import HPL, reduce_known_hpl


@lru_cache(maxsize=4096)
def _extract(word, z, leading, trailing):
    if not word:
        return s.S.One
    if len(set(word)) == 1:
        return HPL((word[0],), z)**len(word)/s.factorial(len(word))
    if trailing and word[-1] == 0:
        prefix = word[:-1]
        terms = shuffle_words(prefix, (0,))
        count = terms.pop(word)
        return (_extract(prefix, z, leading, trailing)*HPL((0,), z)
                - sum(n*_extract(w, z, leading, trailing) for w,n in terms.items()))/count
    if leading and word[0] == 1:
        tail = word[1:]
        terms = shuffle_words((1,), tail)
        count = terms.pop(word)
        return (HPL((1,), z)*_extract(tail, z, leading, trailing)
                - sum(n*_extract(w, z, leading, trailing) for w,n in terms.items()))/count
    return HPL(word, z)


def extract_hpl_logs(expression, *, leading_ones=True, trailing_zeros=True):
    """Remove leading 1 / trailing 0 letters by exact shuffle identities.

    Products of log(z) and -log(1-z) remain explicit; other HPL words retain
    their original arguments. No complex branch or endpoint limit is chosen.
    This is not a minimal-basis reduction or an endpoint evaluator.
    """
    if not isinstance(expression, s.Expr):
        raise TypeError("expression must be a SymPy expression")
    if type(leading_ones) is not bool or type(trailing_zeros) is not bool:
        raise TypeError("Log extraction flags must be booleans")
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("Log extraction requires exact finite input")
    replacements = {}
    for atom in expression.atoms(HPL):
        word = _word(atom.word)
        replacements[atom] = _extract(word, atom.argument, leading_ones, trailing_zeros)
    return reduce_known_hpl(s.expand(expression.xreplace(replacements)))
