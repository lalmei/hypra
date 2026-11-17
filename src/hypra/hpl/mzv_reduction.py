"""Small-weight MZV reductions generated from exact double-shuffle relations.

No numerical fitting or upstream table is used. Admissible shuffle and stuffle
products have the same convergent product on the right. For the otherwise
divergent factor zeta(1), subtract the two products before taking the limit:
their common divergent term cancels. Duality and the beta-integral height-one
identity complete the system through weight eight.
"""

from collections import defaultdict
from collections.abc import Mapping
from fractions import Fraction
from functools import lru_cache
from types import MappingProxyType

import sympy as s

from ..symbolic import MZV
from .algebra import shuffle_words
from .stuffle import _stuffle


@lru_cache(maxsize=9)
def _compositions(weight: int) -> tuple[tuple[int, ...], ...]:
    if not weight:
        return ((),)
    return tuple(
        (first,) + tail for first in range(1, weight + 1) for tail in _compositions(weight - first)
    )


def _admissible(weight):
    return tuple(indices for indices in _compositions(weight) if indices[0] >= 2)


def _word(indices):
    return tuple(letter for index in indices for letter in (0,) * (index - 1) + (1,))


def _indices(word):
    result, zeros = [], 0
    for letter in word:
        if letter == 0:
            zeros += 1
        else:
            result.append(zeros + 1)
            zeros = 0
    return tuple(result)


def _shuffle(left, right):
    return {
        _indices(word): count for word, count in shuffle_words(_word(left), _word(right)).items()
    }


def _quasi_shuffle(left, right):
    return {
        tuple(index for index, _ in word): count
        for word, count in _stuffle(
            tuple((index, 1) for index in left), tuple((index, 1) for index in right)
        )
    }


def _pivot_key(indices):
    # Keep zeta(5,3), the one new weight-eight primitive, when it is free.
    # This is a chosen spanning set, not an algebraic-independence assertion.
    return indices == (5, 3), -len(indices), indices


@lru_cache(maxsize=7)
def mzv_weight_reductions(weight: int) -> Mapping[tuple[int, ...], s.Expr]:
    """Generate exact reductions at a requested weight between two and eight."""
    if type(weight) is not int or not 2 <= weight <= 8:
        raise ValueError("Double-shuffle reduction supports weights two through eight")
    # Imported here because the public reducer invokes this bounded solver.
    from .endpoints import _height_one, dual_indices

    def lower(indices):
        return mzv_weight_reductions(sum(indices))[indices]

    rows = []
    for indices in _admissible(weight):
        if len(indices) == 1:
            rows.append(({indices: 1}, s.zeta(weight)))
        elif all(index == 1 for index in indices[1:]):
            rows.append(({indices: 1}, _height_one(indices[0], len(indices))))
        dual = dual_indices(indices)
        if dual != indices:
            rows.append(({indices: 1, dual: -1}, s.S.Zero))
    for split in range(2, weight - 1):
        for left in _admissible(split):
            for right in _admissible(weight - split):
                rhs = lower(left) * lower(right)
                rows.extend(((_shuffle(left, right), rhs), (_quasi_shuffle(left, right), rhs)))
    if weight > 2:
        for indices in _admissible(weight - 1):
            difference = defaultdict(int, _shuffle((1,), indices))
            for word, count in _quasi_shuffle((1,), indices).items():
                difference[word] -= count
            rows.append((dict(difference), s.S.Zero))

    pivots = {}
    for coefficients, rhs in rows:
        row = {indices: Fraction(count) for indices, count in coefficients.items() if count}
        if any(indices[0] < 2 for indices in row):
            raise ArithmeticError("Divergent double-shuffle terms did not cancel")
        while row:
            pivot = min(row, key=_pivot_key)
            scale = row.pop(pivot)
            if pivot not in pivots:
                pivots[pivot] = (
                    {indices: count / scale for indices, count in row.items()},
                    s.expand(rhs / s.Rational(scale)),
                )
                break
            terms, replacement = pivots[pivot]
            rhs -= s.Rational(scale) * replacement
            for indices, count in terms.items():
                remaining = row.get(indices, 0) - scale * count
                if remaining:
                    row[indices] = remaining
                else:
                    row.pop(indices, None)
        else:
            if s.expand(rhs) != 0:
                raise ArithmeticError("Double-shuffle relations produced an inconsistent row")

    @lru_cache(maxsize=64)
    def resolved(indices):
        if indices not in pivots:
            return MZV(indices)
        terms, rhs = pivots[indices]
        return s.expand(
            rhs - s.Add(*(s.Rational(count) * resolved(other) for other, count in terms.items()))
        )

    return MappingProxyType({indices: resolved(indices) for indices in _admissible(weight)})
