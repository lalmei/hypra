"""Exact HPL index and plus/minus kernel notation conversions."""

from itertools import product
from typing import Literal, Self, cast

import sympy as s
from sympy.core.function import ArgumentIndexError
from sympy.core.symbol import Str

from .._typing import IntegerWord, as_integer_tuple
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL
from .algebra import _word

type PMLetter = Literal[0, "plus", "minus"]
type PMWord = tuple[PMLetter, ...] | list[PMLetter] | s.Tuple
_PM_LETTERS: dict[int, PMLetter] = {-1: "minus", 0: 0, 1: "plus"}


def _pm_indices(word: PMWord) -> tuple[int, ...]:
    if not isinstance(word, (tuple, list, s.Tuple)):
        raise TypeError("Plus/minus word must be a sequence")
    indices = []
    for letter in word:
        if isinstance(letter, Str):
            letter = str(letter)
        if (type(letter) is int or isinstance(letter, s.Integer)) and letter == 0:
            indices.append(0)
        elif letter == "plus":
            indices.append(1)
        elif letter == "minus":
            indices.append(-1)
        else:
            raise ValueError("Use zero, 'plus', or 'minus' as kernel letters")
    return tuple(indices)


class HPLPM(s.Function):
    """An HPL whose word retains the plus/minus kernel notation."""

    nargs = 2

    def __new__(cls, word: PMWord, argument: s.Expr, **options: object) -> Self:
        letters = _pm_indices(word)
        encoded = s.Tuple(*(Str(str(_PM_LETTERS[i])) if i else s.S.Zero for i in letters))
        return cast(Self, super().__new__(cls, encoded, argument, **options))

    @classmethod
    def from_indices(cls, indices: tuple[int, ...], argument: s.Expr) -> Self:
        """Construct from indices -1, 0, 1 representing minus, zero, plus."""
        return cls(tuple(_PM_LETTERS[index] for index in indices), argument)

    @property
    def word(self) -> s.Tuple:
        return cast(s.Tuple, self.args[0])

    @property
    def argument(self) -> s.Expr:
        return cast(s.Expr, self.args[1])

    def _eval_rewrite_as_HPL(self, word, argument, **hints):
        return hpl_pm(word, argument)

    def _eval_evalf(self, prec):
        return self.rewrite(HPL)._eval_evalf(prec)

    def _eval_nseries(self, x, n, logx, cdir=0):
        # The word is discrete notation, not a differentiable Function
        # argument. Reuse the ordinary-letter germ and its remainder checks.
        return self.rewrite(HPL)._eval_nseries(x, n, logx, cdir=cdir)

    def _eval_derivative(self, s):
        return self.fdiff() * self.argument.diff(s)

    def fdiff(self, argindex=2):
        if argindex != 2:
            raise ArgumentIndexError(self, argindex)
        indices = _pm_indices(self.word)
        if not indices:
            return s.S.Zero
        z = self.argument
        first, *rest = indices
        kernel = {-1: 2 * z / (1 - z**2), 0: 1 / z, 1: 2 / (1 - z**2)}[first]
        return kernel * (self.from_indices(tuple(rest), z) if rest else s.S.One)


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
    choices = tuple(choices)
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
    return _linear_choices(
        ((0, s.S.One),) if not letter else ((1, s.S.One), (-1, s.Integer(letter)))
        for letter in _pm_indices(word)
    )


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
