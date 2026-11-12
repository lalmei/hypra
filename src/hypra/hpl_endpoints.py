"""Convergent real HPL endpoints at one and exact multiple-zeta identities.

The endpoint is the limit from 0<z<1. No divergent value is assigned. Trailing
zeros are removed by the shuffle identity with H_0(z)=log(z); its product with
any {0,1} HPL tends to zero at z=1 because HPL endpoint growth is at most a
power of log(1-z). Thus the reduction preserves the complete finite limit.

For admissible words (begin 0, end 1), each run 0^(s-1)1 gives an MZV index s.
MZV duality follows from t -> 1-t and reversing the integration simplex.
Weight-four reductions follow by comparing
  zeta(2)^2 = 2 zeta(2,2) + zeta(4)                  (strict stuffle)
and
  zeta(2)^2 = 2 zeta(2,2) + 4 zeta(3,1)              (HPL shuffle),
with zeta(2)^2 = (5/2) zeta(4). These are identities, not fitted tables.
"""

from __future__ import annotations

from functools import lru_cache

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import shuffle_words
from .symbolic import MZV

_MAX_ENDPOINT_WEIGHT = 24


def _word(word):
    if not isinstance(word, (tuple, list, s.Tuple)):
        raise TypeError("Endpoint input must be an HPL word, not an arbitrary expression")
    if any(
        isinstance(letter, bool) or not isinstance(letter, (int, s.Integer)) or letter not in (0, 1)
        for letter in word
    ):
        raise UnsupportedExpressionError("This endpoint evaluator supports only the {0,1} alphabet")
    if len(word) > _MAX_ENDPOINT_WEIGHT:
        raise UnsupportedExpressionError(f"Endpoint weight exceeds {_MAX_ENDPOINT_WEIGHT}")
    return tuple(int(letter) for letter in word)


def _indices_from_word(word):
    if not word or word[0] != 0 or word[-1] != 1:
        raise UnsupportedExpressionError("An MZV endpoint word must begin with 0 and end with 1")
    indices, zeros = [], 0
    for letter in word:
        if letter == 0:
            zeros += 1
        else:
            indices.append(zeros + 1)
            zeros = 0
    return tuple(indices)


def dual_indices(indices):
    """Convergent MZV duality, derived by reversing and complementing its word."""
    value = MZV(indices)
    if value == 1:
        return ()
    if sum(value.indices) > _MAX_ENDPOINT_WEIGHT:
        raise UnsupportedExpressionError(f"MZV duality weight exceeds {_MAX_ENDPOINT_WEIGHT}")
    word = tuple(letter for index in value.indices for letter in (*([0] * (int(index) - 1)), 1))
    return _indices_from_word(tuple(1 - letter for letter in reversed(word)))


@lru_cache(maxsize=256)
def _height_one(first, depth):
    """Coefficient of x^(first-1)y^depth in 1-Gamma(1-x)Gamma(1-y)/Gamma(1-x-y).

    The convergent height-one MZV generating identity follows from the Euler
    beta integral. Sparse bivariate truncation is exact; every logarithm term
    contains both x and y, so only min(first-1, depth) powers can contribute.
    """
    nx, ny = first - 1, depth
    logarithm = {(i, j): -s.binomial(i + j, i) * s.zeta(i + j) / (i + j)
                 for i in range(1, nx + 1) for j in range(1, ny + 1)}
    power = {(0, 0): s.S.One}
    coefficient = s.S.Zero
    for k in range(1, min(nx, ny) + 1):
        product = {}
        for (i, j), left in power.items():
            for (a, b), right in logarithm.items():
                if i + a <= nx and j + b <= ny:
                    key = (i + a, j + b)
                    product[key] = product.get(key, s.S.Zero) + left * right
        power = product
        coefficient -= power.get((nx, ny), s.S.Zero) / s.factorial(k)
    return s.expand(coefficient)


def reduce_mzv(expression):
    """Apply proved duality, height-one and weight-four identities, opt-in.

    Height-one beta-generating reductions are bounded at weight twelve.
    General irreducible MZVs remain exact symbolic objects. Duality selects the
    representative with smaller depth (then lexicographic order); it does not
    make unproved basis-completeness or transcendence claims.
    """
    if type(expression) is int:
        expression = s.Integer(expression)
    if not isinstance(expression, s.Basic):
        raise TypeError("expression must be a SymPy expression")

    def reduce_one(atom):
        indices = tuple(int(index) for index in atom.indices)
        if len(indices) == 1:
            return s.zeta(indices[0])
        if sum(indices) > _MAX_ENDPOINT_WEIGHT:
            return atom
        dual = dual_indices(indices)
        indices = min((indices, dual), key=lambda word: (len(word), word))
        if len(indices) == 1:
            return s.zeta(indices[0])
        if sum(indices) <= 12 and all(index == 1 for index in indices[1:]):
            return _height_one(indices[0], len(indices))
        if indices == (2, 2):
            return 3 * s.zeta(4) / 4
        return MZV(indices)

    return expression.replace(lambda node: isinstance(node, MZV), reduce_one)


@lru_cache(maxsize=2048)
def _endpoint(word):
    if not word:
        return s.S.One
    if all(letter == 0 for letter in word):
        return s.S.Zero
    # For a leading 1, any further 1 leaves a nonzero (or divergent) tail at
    # one and produces a logarithmic divergence. A pure-zero tail instead
    # vanishes and is integrable: H(1,0;1)=-zeta(2), for example.
    if word[0] == 1 and (len(word) == 1 or 1 in word[1:]):
        raise UnsupportedExpressionError(
            "Divergent HPL endpoint at one; no regularized value is assigned"
        )
    if word[-1] == 1:
        return MZV(_indices_from_word(word))
    prefix = word[:-1]
    shuffled = shuffle_words(prefix, (0,))
    multiplicity = shuffled.pop(word)
    result = -sum((count * _endpoint(other) for other, count in shuffled.items()), s.S.Zero)
    return s.expand(result / multiplicity)


def hpl_at_one(word, *, reduced=True):
    """Evaluate a convergent {0,1} HPL word at z=1 from the real interval.

    The input is the discrete word, e.g. ``hpl_at_one((0,1)) == pi**2/6``.
    Use reduced=False to retain MZV objects without applying proven reductions.
    Empty words give one, nonempty all-zero words give zero. Other alphabets
    and genuinely divergent endpoints raise UnsupportedExpressionError.
    This is not a limit engine for arbitrary rational-HPL combinations.
    """
    if type(reduced) is not bool:
        raise TypeError("reduced must be a Boolean")
    result = _endpoint(_word(word))
    return reduce_mzv(result) if reduced else result
