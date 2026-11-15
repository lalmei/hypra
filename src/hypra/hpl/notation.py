"""Exact HPL index and plus/minus kernel notation conversions."""

from itertools import product
from typing import Literal

import sympy as s

from .._typing import IntegerWord, as_integer_tuple
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL
from .algebra import _word

type PMLetter = Literal[0, "plus", "minus"]
type PMWord = tuple[PMLetter, ...] | list[PMLetter] | s.Tuple


def expand_indices(indices: IntegerWord) -> tuple[int, ...]:
    """Expand signed compressed indices.

    Args:
        indices: Each nonzero index becomes `abs(index) - 1` zero letters
            followed by its sign. A standalone zero stays a zero.

    Returns:
        Expanded letters in their original order.
    """
    if not isinstance(indices, (tuple, list, s.Tuple)) or any(
        isinstance(n, bool) or not isinstance(n, (int, s.Integer)) for n in indices
    ):
        raise TypeError("Compressed indices must be exact integers")
    indices = as_integer_tuple(indices)
    if sum(max(1, abs(n)) for n in indices) > 24:
        raise UnsupportedExpressionError("Index conversion weight exceeds24")
    return tuple(
        letter
        for n in indices
        for letter in ((0,) if n == 0 else (0,) * (abs(n) - 1) + (int(s.sign(n)),))
    )


def compress_indices(word: IntegerWord) -> tuple[int, ...]:
    """Compress an expanded word into signed indices.

    Args:
        word: Expanded letters. Zero runs attach to the next nonzero letter.

    Returns:
        Compressed indices with any final zero letters retained.
    """
    word = _word(word)
    result, zeros = [], 0
    for letter in word:
        if letter == 0:
            zeros += 1
        else:
            result.append(letter * (zeros + 1))
            zeros = 0
    return tuple(result) + (0,) * zeros


def _linear_choices(choices):
    if len(choices) > 12:
        raise UnsupportedExpressionError("Plus/minus conversion weight exceeds12")
    result = {}
    for terms in product(*choices):
        word = tuple(letter for letter, _ in terms)
        result[word] = s.prod(coefficient for _, coefficient in terms)
    return result


def pm_to_letters(word: PMWord) -> dict[tuple[int, ...], s.Expr]:
    """Expand the kernels `f_plus=f_1+f_-1` and `f_minus=f_1-f_-1`.

    Args:
        word: Sequence containing zero, `"plus"`, or `"minus"`.

    Returns:
        Coefficients indexed by expanded HPL words.
    """
    if not isinstance(word, (tuple, list, s.Tuple)):
        raise TypeError("Plus/minus word must be a sequence")
    choices = []
    for letter in word:
        if type(letter) is int and letter == 0 or isinstance(letter, s.Integer) and letter == 0:
            choices.append(((0, s.S.One),))
        elif letter == "plus":
            choices.append(((1, s.S.One), (-1, s.S.One)))
        elif letter == "minus":
            choices.append(((1, s.S.One), (-1, -s.S.One)))
        else:
            raise ValueError("Use zero, 'plus', or 'minus' as kernel letters")
    return _linear_choices(choices)


def letters_to_pm(word: IntegerWord) -> dict[tuple[PMLetter, ...], s.Expr]:
    """Expand the kernels `f_±1=(f_plus ± f_minus)/2`.

    Args:
        word: Expanded HPL letters.

    Returns:
        Coefficients indexed by plus/minus words.
    """
    choices = []
    for letter in _word(word):
        choices.append(
            ((0, s.S.One),)
            if letter == 0
            else (("plus", s.S.Half), ("minus", s.Rational(letter, 2)))
        )
    return _linear_choices(choices)


def hpl_pm(word: PMWord, argument: s.Expr) -> s.Expr:
    """Construct plus/minus notation as an exact sum of ordinary HPL objects.

    Args:
        word: Sequence containing zero, `"plus"`, or `"minus"`.
        argument: SymPy argument of each HPL in the sum.
    """
    if not isinstance(argument, s.Expr):
        raise TypeError("Argument must be a SymPy expression")
    return sum(
        (
            coefficient * (HPL(letters, argument) if letters else s.S.One)
            for letters, coefficient in pm_to_letters(word).items()
        ),
        s.S.Zero,
    )
