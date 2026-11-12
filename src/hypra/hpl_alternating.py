"""Exact convergent alternating HPL endpoints, without regularization tables."""
from functools import lru_cache

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import shuffle_words
from .hpl_endpoints import reduce_mzv
from .hpl_numeric import _word
from .symbolic import MZV, EulerSum


def reduce_euler_sum(expression):
    """Proved eta, MZV, repeated-letter and finite-cutoff Newton identities."""
    def one(atom):
        indices, colors = tuple(atom.indices), tuple(atom.colors)
        if len(indices) > 1 and len(set(indices)) == len(set(colors)) == 1:
            # At finite cutoff these are elementary symmetric polynomials in
            # x_n=c**n/n**s. Newton identities hold before any infinite limit.
            # Every power sum converges for our valid inputs, even s=1,c=-1.
            index, color = indices[0], colors[0]
            elementary = [s.S.One]
            powers = [None] + [one(EulerSum((k * index,), (color ** k,)))
                               for k in range(1, len(indices) + 1)]
            for depth in range(1, len(indices) + 1):
                elementary.append(s.expand(sum(
                    (-1) ** (k - 1) * elementary[depth - k] * powers[k]
                    for k in range(1, depth + 1)
                ) / depth))
            return elementary[-1]
        if all(c == 1 for c in colors):
            return reduce_mzv(MZV(indices))
        if len(indices) == 1:
            n = indices[0]
            return -s.log(2) if n == 1 else -(1 - 2 ** (1 - n)) * s.zeta(n)
        # H((-1)^d;1) = log(2)^d/d!, while its series carries (-1)^d.
        if all(n == 1 for n in indices) and colors == (-1,) + (1,) * (len(indices) - 1):
            return (-1) ** len(indices) * s.log(2) ** len(indices) / s.factorial(len(indices))
        return atom
    if type(expression) is int:
        expression = s.Integer(expression)
    if not isinstance(expression, s.Basic):
        raise TypeError("expression must be a SymPy expression, not source text")
    return expression.replace(lambda node: isinstance(node, EulerSum), one)


@lru_cache(maxsize=2048)
def _endpoint(word):
    if not word:
        return s.S.One
    if all(letter == 0 for letter in word):
        return s.S.Zero
    # Kernels have fixed positive sign on (0,1); a nonzero-letter tail
    # has a nonzero endpoint. Only a pure-zero tail cancels a leading 1 pole.
    if word[0] == 1 and (len(word) == 1 or any(word[1:])):
        raise UnsupportedExpressionError("Divergent HPL endpoint; no regularization")
    if word[-1] == 0:
        shuffled = shuffle_words(word[:-1], (0,))
        multiplicity = shuffled.pop(word)
        return s.expand(-sum((count * _endpoint(w) for w, count in shuffled.items()), s.S.Zero)
                        / multiplicity)
    indices, colors, zeros, previous, factor = [], [], 0, 1, 1
    for letter in word:
        if letter == 0:
            zeros += 1
        else:
            indices.append(zeros + 1)
            colors.append(previous * letter)
            previous = letter
            factor *= letter
            zeros = 0
    return factor * EulerSum(indices, colors)


def hpl_alternating_at_one(word, *, reduced=True):
    """Limit from real 0<z<1 for expanded 0,+1,-1 words, weight <=24.

    Conditional sums use the outer-index cutoff, not arbitrary reordered sums.
    Mixed trailing zeros are removed by shuffle with log(z), which vanishes
    at one even multiplied by any finite power of log(1-z).
    """
    if type(reduced) is not bool:
        raise TypeError("reduced must be Boolean")
    result = _endpoint(_word(word))
    return reduce_euler_sum(result) if reduced else result
