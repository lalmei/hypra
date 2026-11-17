"""Euler-sum reductions from colored double shuffle and distribution.

The four sign choices project onto even summation indices, giving
sum_{colors} ES(a,b;colors) = 2**(2-a-b) zeta(a,b). Combined with shuffle,
stuffle and their convergent single-zeta(1) difference, this determines every
admissible depth-two colored sum at weights three, five and seven. At all
depths through weight six, the pullback u=t/(1+t) supplies Li_n(1/2), closing
the system in zeta, log(2), and Li_n(1/2), with one additional HPLs6 primitive.
"""

from collections import defaultdict
from collections.abc import Mapping
from fractions import Fraction
from functools import lru_cache
from itertools import product
from types import MappingProxyType

import sympy as s

from ..symbolic import MZV, EulerSum
from .algebra import shuffle_words
from .endpoints import reduce_mzv
from .stuffle import _stuffle


def _hpl_word(indices):
    word, previous, sign = [], 1, 1
    for index, color in indices:
        letter = previous * color
        word.extend([0] * (index - 1) + [letter])
        previous = letter
        sign *= letter
    return tuple(word), sign


def _euler_word(word):
    indices, zeros, previous, sign = [], 0, 1, 1
    for letter in word:
        if letter == 0:
            zeros += 1
        else:
            indices.append((zeros + 1, previous * letter))
            zeros, previous = 0, letter
            sign *= letter
    return tuple(indices), sign


def _colored_shuffle(left, right):
    a, sign_a = _hpl_word(left)
    b, sign_b = _hpl_word(right)
    result = defaultdict(int)
    for word, count in shuffle_words(a, b).items():
        indices, sign = _euler_word(word)
        result[indices] += count * sign * sign_a * sign_b
    return dict(result)


def _depth_one(word):
    index, color = word[0]
    if color == 1:
        return s.zeta(index)
    return -s.log(2) if index == 1 else -(1 - s.Rational(1, 2) ** (index - 1)) * s.zeta(index)


@lru_cache(maxsize=3)
def euler_depth_two_reductions(weight: int) -> Mapping[tuple[tuple[int, int], ...], s.Expr]:
    """Return derived reductions for all convergent depth-two sums at weight 3,5,7."""
    if weight not in (3, 5, 7):
        raise ValueError("Depth-two Euler reduction supports odd weights three through seven")
    words = tuple(
        ((split, first), (weight - split, second))
        for split in range(1, weight)
        for first, second in product((-1, 1), repeat=2)
        if (split, first) != (1, 1)
    )
    rows = [
        ({word: 1}, reduce_mzv(MZV(tuple(index for index, _ in word))))
        for word in words
        if all(color == 1 for _, color in word)
    ]
    for split in range(1, weight):
        for first, second in product((-1, 1), repeat=2):
            left, right = ((split, first),), ((weight - split, second),)
            if left == ((1, 1),) or right == ((1, 1),):
                continue
            value = _depth_one(left) * _depth_one(right)
            rows.extend(
                ((_colored_shuffle(left, right), value), (dict(_stuffle(left, right)), value))
            )
    for color in (-1, 1):
        right = ((weight - 1, color),)
        difference = defaultdict(int, _colored_shuffle(((1, 1),), right))
        for word, count in _stuffle(((1, 1),), right):
            difference[word] -= count
        rows.append((dict(difference), s.S.Zero))
    for split in range(2, weight):
        rows.append(
            (
                {
                    ((split, first), (weight - split, second)): 1
                    for first, second in product((-1, 1), repeat=2)
                },
                s.Rational(2) ** (2 - weight) * reduce_mzv(MZV((split, weight - split))),
            )
        )
    return MappingProxyType(_solve(rows, words, depth_one_known=True))


def _scalar_polynomial(expression):
    result = {}
    for term in s.Add.make_args(s.expand(expression)):
        coefficient, monomial = term.as_coeff_Mul(rational=True)
        if coefficient:
            if not isinstance(coefficient, s.Rational):
                raise ArithmeticError("Expected rational coefficients in the constant basis")
            result[monomial] = Fraction(int(coefficient.p), int(coefficient.q))
    return result


def _subtract(target, scale, source):
    for term, coefficient in source.items():
        remaining = target.get(term, 0) - scale * coefficient
        if remaining:
            target[term] = remaining
        else:
            target.pop(term, None)


def _solve(rows, words, *, depth_one_known=False):
    pivots = {}
    for coefficients, rhs in rows:
        row = {word: Fraction(count) for word, count in coefficients.items() if count}
        for word in list(row):
            if depth_one_known and len(word) == 1:
                rhs -= s.Rational(row.pop(word)) * _depth_one(word)
        scalar = _scalar_polynomial(rhs)
        if any(word[0] == (1, 1) for word in row):
            raise ArithmeticError("Divergent colored double-shuffle terms did not cancel")
        while row:
            pivot = min(row, key=lambda word: (-len(word), word))
            scale = row.pop(pivot)
            if pivot not in pivots:
                pivots[pivot] = (
                    {word: count / scale for word, count in row.items()},
                    {term: coefficient / scale for term, coefficient in scalar.items()},
                )
                break
            terms, replacement = pivots[pivot]
            _subtract(scalar, scale, replacement)
            for word, count in terms.items():
                remaining = row.get(word, 0) - scale * count
                if remaining:
                    row[word] = remaining
                else:
                    row.pop(word, None)
        else:
            if scalar:
                raise ArithmeticError(
                    "Colored double-shuffle relations produced an inconsistent row"
                )
    if len(pivots) != len(words):
        raise ArithmeticError("Colored double-shuffle system did not span the requested weight")

    @lru_cache(maxsize=324)
    def resolved(word):
        terms, rhs = pivots[word]
        result = rhs.copy()
        for other, count in terms.items():
            _subtract(result, count, resolved(other))
        return result

    return {
        word: s.Add(*(s.Rational(count) * term for term, count in resolved(word).items()))
        for word in words
    }


@lru_cache(maxsize=7)
def _colored_compositions(weight: int) -> tuple[tuple[tuple[int, int], ...], ...]:
    if not weight:
        return ((),)
    return tuple(
        ((first, color),) + tail
        for first in range(1, weight + 1)
        for color in (-1, 1)
        for tail in _colored_compositions(weight - first)
    )


def _admissible(weight):
    return tuple(word for word in _colored_compositions(weight) if word[0] != (1, 1))


@lru_cache(maxsize=6)
def euler_weight_reductions(weight: int) -> Mapping[tuple[tuple[int, int], ...], s.Expr]:
    """Generate every convergent colored Euler reduction through weight six.

    At weight six, ES((5,1),(-1,-1))=HPLs6-zeta(6) is retained as the
    additional primitive. This chooses a spanning basis without asserting
    algebraic independence of its constants.
    """
    if type(weight) is not int or not 1 <= weight <= 6:
        raise ValueError("General Euler reduction supports weights one through six")
    if weight == 1:
        return MappingProxyType({((1, -1),): -s.log(2)})
    words = _admissible(weight)
    rows = []
    for word in words:
        if len(word) == 1:
            rows.append(({word: 1}, _depth_one(word)))
        elif all(color == 1 for _, color in word):
            rows.append(({word: 1}, reduce_mzv(MZV(tuple(index for index, _ in word)))))
        elif len(word) == 2 and weight % 2:
            rows.append(({word: 1}, euler_depth_two_reductions(weight)[word]))
    for split in range(1, weight):
        for left in _admissible(split):
            for right in _admissible(weight - split):
                rhs = (
                    euler_weight_reductions(split)[left]
                    * euler_weight_reductions(weight - split)[right]
                )
                rows.extend(
                    ((_colored_shuffle(left, right), rhs), (dict(_stuffle(left, right)), rhs))
                )
    for word in _admissible(weight - 1):
        difference = defaultdict(int, _colored_shuffle(((1, 1),), word))
        for indices, count in _stuffle(((1, 1),), word):
            difference[indices] -= count
        rows.append((dict(difference), s.S.Zero))
    for indices in sorted({tuple(index for index, _ in word) for word in words if word[0][0] >= 2}):
        rows.append(
            (
                {
                    tuple(zip(indices, colors, strict=True)): 1
                    for colors in product((-1, 1), repeat=len(indices))
                },
                s.Rational(2) ** (len(indices) - weight) * reduce_mzv(MZV(indices)),
            )
        )

    # Under u=t/(1+t), omega_0(u)=omega_0(t)-omega_-1(t) and
    # omega_1(u)=omega_-1(t). Every term ends in -1, so all endpoint sums
    # converge. This supplies the half-argument periods at weights four to six.
    half_value = defaultdict(int)
    for prefix in product((0, -1), repeat=weight - 1):
        indices, sign = _euler_word(prefix + (-1,))
        half_value[indices] += (-1) ** prefix.count(-1) * sign
    polylog = s.polylog(weight, s.S.Half)
    if weight == 3:
        polylog = 7 * s.zeta(3) / 8 - s.pi**2 * s.log(2) / 12 + s.log(2) ** 3 / 6
    rows.append((dict(half_value), polylog))
    if weight == 6:
        # Choose the one additional period in the shipped weight-six basis.
        # Keeping its exact strict sum avoids introducing a new opaque symbol.
        primitive = ((5, -1), (1, -1))
        rows.append(({primitive: 1}, EulerSum((5, 1), (-1, -1))))
    return MappingProxyType(_solve(rows, words))
