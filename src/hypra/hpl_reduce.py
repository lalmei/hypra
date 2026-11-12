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


def lyndon_factors(word):
    """Chen--Fox--Lyndon factorization, with alphabet order -1 < 0 < 1."""
    word = _word(word)
    factors, i = [], 0
    while i < len(word):
        j, k = i+1, i
        while j < len(word) and word[k] <= word[j]:
            k = i if word[k] < word[j] else k+1
            j += 1
        length = j-k
        while i <= k:
            factors.append(word[i:i+length])
            i += length
    return tuple(factors)


@lru_cache(maxsize=4096)
def _lyndon(word, z):
    if not word:
        return s.S.One
    factors = lyndon_factors(word)
    if len(factors) == 1:
        return HPL(word, z)
    counts = {(): 1}
    for factor in factors:
        merged = {}
        for left, multiplicity in counts.items():
            for right, coefficient in shuffle_words(left, factor).items():
                merged[right] = merged.get(right, 0) + multiplicity*coefficient
        counts = merged
        if len(counts) > 16384:
            raise UnsupportedExpressionError("Lyndon reduction exceeds the word budget")
    diagonal = counts.pop(word)
    if any(other >= word for other in counts):
        raise UnsupportedExpressionError("Lyndon shuffle triangularity check failed")
    product = s.prod(HPL(factor, z) for factor in factors)
    return s.expand((product-sum(count*_lyndon(other,z) for other,count in counts.items()))/diagonal)


def lyndon_reduce(expression):
    """Reduce HPL words to a shuffle-algebra Lyndon basis through weight eight.

    This uses exact triangular shuffle factorization, not precomputed tables.
    It is a shuffle basis only, not a claim of algebraic independence after
    argument transformations or evaluation at special points.
    """
    if not isinstance(expression, s.Expr):
        raise TypeError("expression must be a SymPy expression")
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("Lyndon reduction requires exact finite input")
    replacements = {}
    for atom in expression.atoms(HPL):
        word = _word(atom.word)
        if len(word) > 8:
            raise UnsupportedExpressionError("Lyndon reduction currently supports weight at most eight")
        replacements[atom] = _lyndon(word, atom.argument)
    return s.expand(expression.xreplace(replacements))
