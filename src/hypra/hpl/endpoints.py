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
from typing import cast

import sympy as s

from .._typing import ExactScalar, IntegerWord, as_integer_tuple
from ..errors import UnsupportedExpressionError
from ..symbolic import MZV
from .algebra import shuffle_words

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
    return as_integer_tuple(word)


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
    if not isinstance(value, MZV):
        return ()
    indices = as_integer_tuple(value.indices)
    if sum(indices) > _MAX_ENDPOINT_WEIGHT:
        raise UnsupportedExpressionError(f"MZV duality weight exceeds {_MAX_ENDPOINT_WEIGHT}")
    word = tuple(letter for index in indices for letter in (*([0] * (index - 1)), 1))
    return _indices_from_word(tuple(1 - letter for letter in reversed(word)))


@lru_cache(maxsize=256)
def _height_one(first, depth):
    """Coefficient of x^(first-1)y^depth in 1-Gamma(1-x)Gamma(1-y)/Gamma(1-x-y).

    The convergent height-one MZV generating identity follows from the Euler
    beta integral. Sparse bivariate truncation is exact; every logarithm term
    contains both x and y, so only min(first-1, depth) powers can contribute.
    """
    nx, ny = first - 1, depth
    logarithm = {
        (i, j): -s.binomial(i + j, i) * s.zeta(i + j) / (i + j)
        for i in range(1, nx + 1)
        for j in range(1, ny + 1)
    }
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


@lru_cache(maxsize=64)
def _repeated_index(index, depth):
    """Newton identities for elementary sums in x_n=1/n**index."""
    elementary: list[s.Expr] = [s.S.One]
    for degree in range(1, depth + 1):
        elementary.append(
            s.expand(
                sum(
                    (-1) ** (k - 1) * elementary[degree - k] * s.zeta(k * index)
                    for k in range(1, degree + 1)
                )
                / degree
            )
        )
    return elementary[-1]


def reduce_mzv(expression: ExactScalar) -> s.Expr:
    """Apply proved duality, height-one and bounded double-shuffle identities.

    Double-shuffle reduction spans all convergent MZVs through weight eight
    using ordinary zeta products and zeta(5,3). Height-one beta-generating
    reductions extend through weight twelve. Repeated indices and Euler's
    zeta(n,1), Zagier's (3,1)^n and the prefixed (2,(1,3)^n) identities extend
    through the endpoint weight limit of 24.
    Higher MZVs remain exact objects.
    Duality selects the
    representative with smaller depth (then lexicographic order); it does not
    make unproved basis-completeness or transcendence claims.

    Args:
        expression: Exact symbolic expression containing multiple zeta values.
    """
    if type(expression) is int:
        expression = s.Integer(expression)
    if not isinstance(expression, s.Basic):
        raise TypeError("expression must be a SymPy expression")

    def reduce_one(atom):
        indices = as_integer_tuple(atom.indices)
        if len(indices) == 1:
            return s.zeta(indices[0])
        if sum(indices) > _MAX_ENDPOINT_WEIGHT:
            return atom
        dual = dual_indices(indices)
        indices = min((indices, dual), key=lambda word: (len(word), word))
        if len(indices) == 1:
            return s.zeta(indices[0])
        # Borwein--Bradley--Broadhurst--Lisonek, arXiv:math/9910045,
        # Example 2.2 / Section 11 proves the Zagier family by a product of
        # Gauss hypergeometric generating functions.
        if len(indices) % 2 == 0 and indices == (3, 1) * (len(indices) // 2):
            repeats = len(indices) // 2
            return 2 * s.pi ** (4 * repeats) / s.factorial(4 * repeats + 2)
        # Example 2.3 (also 6.4) of the same paper records the proved
        # Bowman--Bradley prefixed family. Repeated-four factors below are
        # generated by Newton identities, independently of any source table.
        if len(indices) % 2 and indices == (2,) + (1, 3) * (len(indices) // 2):
            repeats = len(indices) // 2
            return s.expand(
                s.Rational(1, 4) ** repeats
                * sum(
                    (-1) ** k
                    * _repeated_index(4, repeats - k)
                    * (
                        (4 * k + 1) * s.zeta(4 * k + 2)
                        - 4
                        * sum(
                            s.zeta(4 * j - 1) * s.zeta(4 * k - 4 * j + 3) for j in range(1, k + 1)
                        )
                    )
                    for k in range(repeats + 1)
                )
            )
        if indices and len(set(indices)) == 1:
            return _repeated_index(indices[0], len(indices))
        if len(indices) == 2 and indices[1] == 1:
            first = indices[0]
            return s.expand(
                first * s.zeta(first + 1) / 2
                - sum(s.zeta(first - k) * s.zeta(k + 1) for k in range(1, first - 1)) / 2
            )
        if indices and sum(indices) <= 12 and all(index == 1 for index in indices[1:]):
            return _height_one(indices[0], len(indices))
        if sum(indices) <= 8:
            from .mzv_reduction import mzv_weight_reductions

            return mzv_weight_reductions(sum(indices))[indices]
        return MZV(indices)

    return cast(s.Expr, expression.replace(lambda node: isinstance(node, MZV), reduce_one))


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


def hpl_at_one(word: IntegerWord, *, reduced: bool = True) -> s.Expr:
    """Evaluate a convergent {0,1} HPL word at z=1 from the real interval.

    The input is the discrete word, e.g. ``hpl_at_one((0,1)) == pi**2/6``.
    Use reduced=False to retain MZV objects without applying proven reductions.
    Empty words give one, nonempty all-zero words give zero. Other alphabets
    and genuinely divergent endpoints raise UnsupportedExpressionError.
    This is not a limit engine for arbitrary rational-HPL combinations.

    Args:
        word: Expanded binary word with a convergent endpoint.
        reduced: Apply supported exact MZV reductions.
    """
    if type(reduced) is not bool:
        raise TypeError("reduced must be a Boolean")
    result = _endpoint(_word(word))
    return reduce_mzv(result) if reduced else result
