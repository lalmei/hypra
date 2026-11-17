"""At-i constant bases generated from shuffle, Möbius and square identities.

Only the retained primitive word names are chosen to agree with HPL 2.0.
Every coefficient is derived by rational elimination; no formula table or
numerical integer-relation fit is used. The primitive atoms are ordinary
HPLPM(word, I) objects and keep their independent numerical evaluation.
"""

from collections.abc import Mapping
from fractions import Fraction
from functools import lru_cache
from itertools import product
from types import MappingProxyType

import sympy as s

from .algebra import shuffle_words
from .alternating import hpl_alternating_at_one
from .euler_reduction import _scalar_polynomial, _subtract
from .notation import HPLPM, pm_to_letters
from .regularization import hpl_endpoint_polynomial

_ALPHABET = (-1, 0, 1)  # minus, zero, plus
_ENDPOINT_LOG = s.Symbol("at_i_endpoint_log")
_MOBIUS_FORMS = {
    0: ((1, -1),),
    1: ((0, -1),),
    -1: ((0, -1), (1, 1), (-1, -1)),
}

# Identifiers of the free HPLI constants in the pinned HPLatI4/5/6 sources.
# They specify coordinates only; the reduction coefficients are generated.
_PRIMITIVES = {
    4: ((0, -1, 1, -1),),
    5: (
        (-1, 1, -1, -1, -1),
        (-1, 1, -1, -1, 1),
        (-1, 1, -1, 1, 1),
        (-1, 1, 1, -1, -1),
        (-1, 1, 1, 1, -1),
        (0, -1, 1, -1, -1),
        (0, 1, -1, -1, -1),
        (0, 1, 1, -1, -1),
        (0, 1, 1, -1, 1),
    ),
    6: (
        (-1, 1, -1, -1, -1, -1),
        (-1, 1, -1, -1, -1, 1),
        (-1, 1, -1, -1, 1, 1),
        (-1, 1, -1, 1, 1, 1),
        (-1, 1, 1, -1, -1, -1),
        (-1, 1, 1, -1, -1, 1),
        (-1, 1, 1, 1, -1, -1),
        (-1, 1, 1, 1, -1, 1),
        (0, -1, -1, 1, -1, -1),
        (0, -1, 1, -1, -1, -1),
        (0, -1, 1, 1, 1, -1),
        (0, 1, -1, -1, -1, -1),
        (0, 1, -1, -1, 0, -1),
        (0, 1, -1, 0, -1, -1),
        (0, 1, -1, 1, 0, -1),
        (0, 1, 1, -1, -1, -1),
        (0, 1, 1, -1, -1, 1),
        (0, 1, 1, -1, 0, -1),
        (0, 1, 1, -1, 0, 1),
        (0, 1, 1, -1, 1, -1),
        (0, 1, 1, 0, -1, -1),
        (0, 1, 1, 1, -1, -1),
        (0, 1, 1, 1, 0, -1),
    ),
}


@lru_cache(maxsize=1093)
def _endpoint_constant(word):
    # M(z)=(1-z)/(1+z) has -log(1-M(z))=-log(z)-log(2)+o(1).
    # Matching the constant term at the tangential origin fixes the pullback.
    letters = pm_to_letters(HPLPM.from_indices(word, s.I).word)
    return s.expand(
        s.Add(
            *(
                coefficient
                * hpl_endpoint_polynomial(ordinary, _ENDPOINT_LOG).subs(_ENDPOINT_LOG, -s.log(2))
                for ordinary, coefficient in letters.items()
            )
        )
    )


@lru_cache(maxsize=1093)
def _mobius_terms(word):
    if not word:
        return (((), s.S.One),)
    result = {}
    for tail, coefficient in _mobius_terms(word[1:]):
        for letter, factor in _MOBIUS_FORMS[word[0]]:
            key = (letter,) + tail
            result[key] = result.get(key, s.S.Zero) + factor * coefficient
    result[()] = _endpoint_constant(word)
    return tuple((key, s.expand(value)) for key, value in result.items() if value != 0)


def _trailing_zeros(word):
    count = 0
    for letter in reversed(word):
        if letter:
            break
        count += 1
    return count


@lru_cache(maxsize=127)
def _square_value(word):
    # Under u=z², omega_0(u)=2*omega_0(z), omega_1(u)=omega_minus(z).
    # At z=i, u=-1 on its upper boundary. Trailing zeros acquire log(-1)=i*pi.
    count = _trailing_zeros(word)
    return s.expand(
        s.Rational(1, 2) ** word.count(0)
        * (-1) ** word.count(-1)
        * s.Add(
            *(
                (s.I * s.pi) ** k / s.factorial(k) * hpl_alternating_at_one(word[: len(word) - k])
                for k in range(count + 1)
            )
        )
    )


def _fractional_row(weight):
    # Li_n(z/(z-1)) has kernel pullbacks omega_0 -> omega_0+omega_1,
    # omega_1 -> -omega_1. Ordinary omega_1=(omega_plus+omega_minus)/2.
    # Continuing z from (0,1) to i stays off the original function's cut.
    result = {}
    for prefix in product(_ALPHABET, repeat=weight - 1):
        coefficient = -(s.Rational(1, 2) ** (1 + sum(letter != 0 for letter in prefix)))
        for final in (-1, 1):
            result[prefix + (final,)] = coefficient
    return result


@lru_cache(maxsize=7)
def at_i_weight_reductions(weight: int) -> Mapping[tuple[int, ...], s.Expr]:
    """Canonical PM at-i values through weight six, with immutable cached results."""
    if type(weight) is not int or not 0 <= weight <= 6:
        raise ValueError("The at-i basis supports weights zero through six")
    if weight == 0:
        return MappingProxyType({(): s.S.One})

    def lower(word):
        return at_i_weight_reductions(len(word))[word]

    words = tuple(product(_ALPHABET, repeat=weight))
    primitives = frozenset(_PRIMITIVES.get(weight, ()))
    pivots = {}

    def add(coefficients, rhs):
        row = {
            word: Fraction(coefficient) for word, coefficient in coefficients.items() if coefficient
        }
        scalar = _scalar_polynomial(rhs)
        while row:
            pivot = min(
                row,
                key=lambda word: (word in primitives, word[0] != 1 and word[-1] != 0, word),
            )
            scale = row.pop(pivot)
            if pivot not in pivots:
                pivots[pivot] = (
                    {word: coefficient / scale for word, coefficient in row.items()},
                    {term: coefficient / scale for term, coefficient in scalar.items()},
                )
                break
            terms, replacement = pivots[pivot]
            _subtract(scalar, scale, replacement)
            _subtract(row, scale, terms)
        else:
            if scalar:
                raise ArithmeticError("At-i constant identities produced an inconsistent row")

    for split in range(1, weight // 2 + 1):
        for left in product(_ALPHABET, repeat=split):
            for right in product(_ALPHABET, repeat=weight - split):
                add(shuffle_words(left, right), lower(left) * lower(right))

    for word in words:
        parity = (-1) ** word.count(1)
        row, rhs = {word: parity}, s.S.Zero
        for target, coefficient in _mobius_terms(word):
            if len(target) == weight:
                row[target] = row.get(target, 0) - coefficient
            else:
                rhs += coefficient * lower(target)
        # M(i)=-i is reached through the lower quadrant. Negation sends each
        # plus letter to its negative; the tangential log shifts by -i*pi.
        rhs -= parity * s.Add(
            *(
                (-s.I * s.pi) ** k / s.factorial(k) * lower(word[:-k])
                for k in range(1, _trailing_zeros(word) + 1)
            )
        )
        add(row, rhs)

    for word in product((-1, 0), repeat=weight):
        add({word: 1}, _square_value(word))
    from .special import _quarter_period_polylog

    classical = (
        s.I * s.pi / 2
        if weight == 1
        else _quarter_period_polylog(weight, 1) - _quarter_period_polylog(weight, -1)
    )
    add({(0,) * (weight - 1) + (1,): 1}, classical)
    if weight in (3, 4):
        add(_fractional_row(weight), s.polylog(weight, (1 - s.I) / 2))
    if set(words).difference(pivots) != primitives:
        raise ArithmeticError("At-i identities did not span the chosen primitive basis")

    @lru_cache(maxsize=729)
    def resolved(word):
        if word in primitives:
            return {HPLPM.from_indices(word, s.I): Fraction(1)}
        terms, rhs = pivots[word]
        result = rhs.copy()
        for other, coefficient in terms.items():
            _subtract(result, coefficient, resolved(other))
        return result

    return MappingProxyType(
        {
            word: s.Add(*(s.Rational(count) * term for term, count in resolved(word).items()))
            for word in words
        }
    )
