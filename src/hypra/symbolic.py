"""Lossless harmonic-polylogarithm representation.

``HPL((a1, ..., an), z)`` uses the expanded Remiddi--Vermaseren alphabet
``ai in {-1, 0, 1}`` with kernels ``1/(1+t)``, ``1/t`` and ``1/(1-t)``.
In particular H_1(z) = -log(1-z), not log(1-z). The object deliberately
retains the symbolic function, including its word, without branch-changing
simplifications. Its analytic continuation is the principal continuation
from the real interval 0 < z < 1 unless an external calculation specifies
otherwise. It does not encode a separate continuation path.
"""

from __future__ import annotations

from sympy import Basic, Function, Integer, Tuple, factorial, log, polylog
from sympy.core.function import ArgumentIndexError


class HPL(Function):
    """A harmonic polylogarithm with an explicit, expanded letter word.

    Use :meth:`from_compressed` for HypExp/HPL's compressed integer indices:
    ``HPL.from_compressed((2, -1), z) == HPL((0, 1, -1), z)``.
    A compressed zero is one zero letter. No numerical or branch-dependent
    reduction is performed automatically.
    """

    nargs = 2

    @classmethod
    def eval(cls, word, argument):
        if not isinstance(word, Tuple):
            raise TypeError("HPL word must be a tuple of letters -1, 0, or 1")
        if any(not isinstance(letter, Integer) or letter not in (-1, 0, 1) for letter in word):
            raise ValueError("HPL requires expanded letters -1, 0, or 1; use from_compressed")

    @classmethod
    def from_compressed(cls, word, argument):
        """Expand signed integer indices without changing their order."""
        if not isinstance(word, (tuple, list, Tuple)):
            raise TypeError("Compressed HPL indices must be a tuple or list of integers")
        letters = []
        for index in word:
            if isinstance(index, bool) or not isinstance(index, (int, Integer)):
                raise TypeError("Compressed HPL indices must be exact integers")
            index = int(index)
            if abs(index) > 10000:
                raise ValueError("Compressed HPL index exceeds the supported word length")
            if index == 0:
                letters.append(Integer(0))
            else:
                letters.extend([Integer(0)] * (abs(index) - 1))
                letters.append(Integer(1 if index > 0 else -1))
            if len(letters) > 10000:
                raise ValueError("HPL word exceeds the supported length")
        return cls(Tuple(*letters), argument)

    @property
    def word(self):
        """Expanded letters, as a SymPy Tuple."""
        return self.args[0]

    @property
    def argument(self):
        return self.args[1]

    def fdiff(self, argindex=2):
        """Differentiate with respect to the argument, not the discrete word.

        The kernel is 1/z for a leading 0, 1/(1-z) for 1 and
        1/(1+z) for -1. The empty word is the constant function one.
        """
        if argindex != 2:
            raise ArgumentIndexError(self, argindex)
        if not self.word:
            return Integer(0)
        z = self.argument
        letter = self.word[0]
        kernel = {Integer(0): 1 / z, Integer(1): 1 / (1 - z), Integer(-1): 1 / (1 + z)}[
            letter
        ]
        tail = HPL(self.word[1:], z) if len(self.word) > 1 else Integer(1)
        return kernel * tail

    def _eval_derivative(self, symbol):
        # Function's generic implementation tries to differentiate Tuple(word),
        # but the letters are discrete constants. Only the argument varies.
        return self.fdiff(2) * self.argument.diff(symbol)

    def _latex(self, printer):
        word = ",".join(printer._print(letter) for letter in self.word)
        return rf"H_{{{word}}}\!\left({printer._print(self.argument)}\right)"



def reduce_known_hpl(expression):
    """Reduce selected HPL words to exact logarithms and classical polylogs.

    This is deliberately opt-in: construction/serialization retain HPL objects.
    Conventions are the principal branches continued from 0 < z < 1, with the
    same side of cuts in both HPL and output functions. This routine performs
    no argument transformations, no power/log simplification and no assertion
    that the expression is regular at a branch point.

    Supported identities are the empty word (one), all-zero, all-one,
    all-minus-one words, and zero strings followed by one. Mixed words outside
    these patterns remain symbolic. Formulae hold away from singular endpoints
    and extend there only when the corresponding limit exists.
    """
    if type(expression) is int:
        expression = Integer(expression)
    if not isinstance(expression, Basic):
        raise TypeError("expression must be a SymPy expression")

    def reduce_one(atom):
        word, z = atom.word, atom.argument
        weight = len(word)
        if not word:
            return Integer(1)
        if all(letter == 0 for letter in word):
            return log(z) ** weight / factorial(weight)
        if all(letter == 1 for letter in word):
            return (-log(1 - z)) ** weight / factorial(weight)
        if all(letter == -1 for letter in word):
            return log(1 + z) ** weight / factorial(weight)
        if word[-1] == 1 and all(letter == 0 for letter in word[:-1]):
            return polylog(weight, z)
        return atom

    return expression.replace(lambda node: isinstance(node, HPL), reduce_one)
