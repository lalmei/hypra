"""Minimal HPL bases generated from shuffle relations."""

from fractions import Fraction
from functools import lru_cache
from itertools import product

import sympy as s

from ..errors import UnsupportedExpressionError
from ..symbolic import HPL
from .algebra import _word, shuffle_words
from .notation import HPLPM, _pm_indices

_ALPHABET = (-1, 0, 1)
_ARGUMENT = s.Dummy("argument")


def _pivot_key(word):
    # Extract endpoint logarithms first. In the plus/minus alphabet, the same
    # ordering eliminates leading plus and trailing zero.
    return word[0] != 1 and word[-1] != 0, word


@lru_cache(maxsize=5)
def _shuffle_pivots(weight):
    """Solve each homogeneous weight over products of lower-weight words."""
    pivots = {}
    for split in range(1, weight // 2 + 1):
        left_words = product(_ALPHABET, repeat=split)
        right_words = product(_ALPHABET, repeat=weight - split)
        for left, right in product(left_words, right_words):
            row = {word: Fraction(count) for word, count in shuffle_words(left, right).items()}
            rhs = HPL(left, _ARGUMENT) * HPL(right, _ARGUMENT)
            while row:
                pivot = min(row, key=_pivot_key)
                scale = row.pop(pivot)
                if pivot not in pivots:
                    pivots[pivot] = (
                        {word: count / scale for word, count in row.items()},
                        rhs / s.Rational(scale),
                    )
                    break
                terms, replacement = pivots[pivot]
                rhs -= s.Rational(scale) * replacement
                for word, count in terms.items():
                    remaining = row.get(word, 0) - scale * count
                    if remaining:
                        row[word] = remaining
                    else:
                        row.pop(word, None)
    return pivots


@lru_cache(maxsize=1092)
def _minimal_word(word):
    if not word:
        return s.S.One
    if len(word) == 1:
        return HPL(word, _ARGUMENT)
    pivots = _shuffle_pivots(len(word))
    if word not in pivots:
        return HPL(word, _ARGUMENT)
    terms, rhs = pivots[word]
    lower = rhs.xreplace({atom: _minimal_word(_word(atom.word)) for atom in rhs.atoms(HPL)})
    return s.expand(
        lower - sum(s.Rational(count) * _minimal_word(other) for other, count in terms.items())
    )


def reduce_hpl_to_minimal_set(expression: s.Expr) -> s.Expr:
    """Reduce ordinary and plus/minus HPLs to the HPL 2.0 shuffle bases.

    Supports weights zero through six. Use HPLPM to retain plus/minus letters
    in the returned basis. The two alphabets are reduced independently.
    """
    if not isinstance(expression, s.Expr):
        raise TypeError("expression must be a SymPy expression")
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("Minimal-set reduction requires exact finite input")
    replacements = {}
    for atom in expression.atoms(HPL, HPLPM):
        if len(atom.word) > 6:
            raise UnsupportedExpressionError("Minimal-set reduction supports weight at most six")
        plus_minus = isinstance(atom, HPLPM)
        template = _minimal_word(_pm_indices(atom.word) if plus_minus else _word(atom.word))
        constructor = HPLPM.from_indices if plus_minus else HPL
        replacements[atom] = template.xreplace(
            {
                component: constructor(_word(component.word), atom.argument)
                for component in template.atoms(HPL)
            }
        )
    return s.expand(expression.xreplace(replacements))
