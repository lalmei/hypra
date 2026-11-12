"""Native canonical hypergeometric expansion by strict finite harmonic sums.

No external tables are used. For a positive-index word ``a`` we use

    H_N(a1,...,ad) = sum_{N >= n1 > ... > nd >= 1} prod_j n_j**(-aj).

The empty word has value one. Products use the positive-merge quasi-shuffle
(stuffle) algebra. For example H_N(1)**2 = 2 H_N(1,1) + H_N(2).
These finite sums are intermediate algebra, not numerical approximations.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from functools import lru_cache

import sympy as s

from .symbolic import HPL

Word = tuple[int, ...]
HarmonicPolynomial = dict[Word, s.Expr]


def _word(value) -> Word:
    if not isinstance(value, (tuple, list, s.Tuple)):
        raise TypeError("A harmonic word must be a tuple/list of positive integer indices")
    if any(
        isinstance(index, bool) or not isinstance(index, (int, s.Integer)) or index <= 0
        for index in value
    ):
        raise ValueError("Strict harmonic words require positive integer indices")
    return tuple(int(index) for index in value)


def _coefficient(value) -> s.Expr:
    if type(value) is int:
        value = s.Integer(value)
    if not isinstance(value, s.Expr) or value.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise ValueError("Harmonic coefficients must be exact finite SymPy expressions")
    return value


def _clean(poly: Mapping) -> HarmonicPolynomial:
    if not isinstance(poly, Mapping):
        raise TypeError("Harmonic polynomials must map words to exact coefficients")
    result = {}
    for word, coefficient in poly.items():
        word, coefficient = _word(word), _coefficient(coefficient)
        value = s.expand(result.get(word, s.S.Zero) + coefficient)
        if value == 0:
            result.pop(word, None)
        else:
            result[word] = value
    return result


@lru_cache(maxsize=32768)
def _quasi_shuffle(left: Word, right: Word) -> tuple[tuple[Word, int], ...]:
    if not left:
        return ((right, 1),)
    if not right:
        return ((left, 1),)
    terms = {}
    for prefix, rest_left, rest_right in (
        (left[0], left[1:], right),
        (right[0], left, right[1:]),
        (left[0] + right[0], left[1:], right[1:]),
    ):
        for tail, multiplicity in _quasi_shuffle(rest_left, rest_right):
            word = (prefix, *tail)
            terms[word] = terms.get(word, 0) + multiplicity
    return tuple(sorted(terms.items()))


def quasi_shuffle(left, right) -> dict[Word, int]:
    """Multiply two strict finite harmonic words, retaining multiplicities.

    The return value is a fresh dictionary; modifying it never changes cached
    products. Merged indices have a positive sign for strict nested sums.
    """
    return dict(_quasi_shuffle(_word(left), _word(right)))


def harmonic_add(*polynomials: Mapping) -> HarmonicPolynomial:
    """Add finite linear combinations of strict harmonic words."""
    result = {}
    for poly in polynomials:
        for word, coefficient in _clean(poly).items():
            result[word] = result.get(word, s.S.Zero) + coefficient
    return _clean(result)


def harmonic_scale(polynomial: Mapping, coefficient) -> HarmonicPolynomial:
    """Scale a harmonic polynomial by an exact symbolic scalar."""
    scalar = _coefficient(coefficient)
    return _clean({word: scalar * value for word, value in _clean(polynomial).items()})


def harmonic_product(left: Mapping, right: Mapping) -> HarmonicPolynomial:
    """Multiply harmonic polynomials using the strict quasi-shuffle algebra."""
    terms = {}
    for left_word, left_coefficient in _clean(left).items():
        for right_word, right_coefficient in _clean(right).items():
            for word, multiplicity in _quasi_shuffle(left_word, right_word):
                terms[word] = terms.get(word, s.S.Zero) + (
                    multiplicity * left_coefficient * right_coefficient
                )
    return _clean(terms)


def harmonic_exp(log_coefficients: Mapping[int, Mapping], order: int) -> list[HarmonicPolynomial]:
    """Coefficients of exp(sum_{m>=1} eps**m L_m), through ``order``.

    ``L_m`` are harmonic polynomials, and must have strictly positive integer
    degree keys. The recurrence is F_0=1 and
    F_n = (1/n) sum_{m=1}^n m L_m * F_{n-m}, where * is quasi-shuffle.
    Missing L_m are zero. No truncation term is silently inferred from data.
    """
    _order(order)
    if not isinstance(log_coefficients, Mapping):
        raise TypeError("Logarithmic coefficients must be a degree-to-polynomial mapping")
    logarithm = {}
    for degree, poly in log_coefficients.items():
        if type(degree) is not int or degree <= 0:
            raise ValueError("Logarithmic degrees must be positive Python integers")
        logarithm[degree] = _clean(poly)
    coefficients = [{(): s.S.One}]
    for degree in range(1, order + 1):
        total = {}
        for k in range(1, degree + 1):
            total = harmonic_add(
                total,
                harmonic_scale(harmonic_product(logarithm.get(k, {}), coefficients[degree - k]), k),
            )
        coefficients.append(harmonic_scale(total, s.Rational(1, degree)))
    return coefficients


def harmonic_word_to_hpl(word, argument):
    """Li_(s1,...,sd)(z) as an expanded 0/1 HPL word, empty word = 1."""
    indices = _word(word)
    if not indices:
        return s.S.One
    letters = tuple(letter for index in indices for letter in (*([0] * (index - 1)), 1))
    return HPL(letters, argument)


def generating_sum(power: int, word, argument):
    """Exact generating function sum_{n>=1} z^n/n^power H_(n-1)(word).

    Initially convergent for |z|<1. The native caller enforces its declared real
    interval 0<z<1 and does not use this identity for a branch-changing rewrite.
    For power=0, summing the outer geometric series produces z/(1-z).
    """
    if type(power) is not int or power < 0:
        raise ValueError("Generating-sum power must be a nonnegative integer")
    word = _word(word)
    if power:
        return harmonic_word_to_hpl((power, *word), argument)
    return argument / (1 - argument) * harmonic_word_to_hpl(word, argument)


def _order(order):
    if type(order) is not int or order < 0:
        raise ValueError("Expansion order must be a nonnegative Python integer")


def _rational(value):
    if type(value) is int:
        value = s.Integer(value)
    if not isinstance(value, s.Rational):
        raise TypeError("Canonical parameter bases must be exact rational numbers")
    return value


def _slope(value, regulator=None):
    """An exact internal polynomial slope, optionally in formal symbols.

    Public input planning remains rational-affine. Internal contiguous
    regularization may introduce formal parameters; coefficients must be exact
    algebraic numbers and dependence on those parameters must be polynomial.
    """
    value = _coefficient(value)
    if regulator is not None and value.has(regulator):
        raise ValueError("A canonical slope must be independent of the expansion regulator")
    symbols = tuple(sorted(value.free_symbols, key=s.default_sort_key))
    if any(type(symbol) is not s.Symbol or symbol.is_commutative is False for symbol in symbols):
        raise ValueError("Formal slopes require ordinary commutative Symbols")
    if not symbols:
        if value.is_algebraic is not True:
            raise ValueError("Constant slopes must be exact algebraic numbers")
        return value
    try:
        polynomial = s.Poly(value, *symbols, extension=True)
    except (s.PolynomialError, s.CoercionFailed) as exc:
        raise ValueError("Internal symbolic slopes must be algebraic polynomials") from exc
    if any(coefficient.is_algebraic is not True for coefficient in polynomial.coeffs()):
        raise ValueError("Internal slope polynomial coefficients must be exact algebraic numbers")
    return s.expand(value)


def canonical_hyper_series(upper_pairs, lower_slopes, z, eps, order):
    """Expand canonical pF(p-1) with upper bases 0/1 and all lower bases 1.

    ``upper_pairs`` contains (base, slope) for base+slope*eps; lower parameters
    are 1+slope*eps. Slopes may be exact algebraic polynomials in internal
    formal parameters (the public planner admits rational slopes only). The result is an exact finite
    polynomial through eps**order, with HPL coefficients. It contains no O term:
    the caller owns the final truncation/provenance and prefactor depth budget.

    For r zero-base uppers and A their slope product, the n>=1 coefficient is

      A eps**r / n**r * exp(sum_{m>=1} (-1)**(m+1) eps**m/m *
        [(sum_upper slopes**m - sum_lower slopes**m)*H_(n-1)(m)
         +(sum_base1 slopes**m - sum_lower slopes**m)/n**m]).

    The harmonic and rational n parts are expanded by formal exponential
    recurrences and then summed by positive-index multiple polylogarithms.
    This algorithm is not an imported coefficient table.
    """
    _order(order)
    if not isinstance(upper_pairs, Sequence) or not isinstance(lower_slopes, Sequence):
        raise TypeError("Canonical parameters must be finite sequences")
    if len(upper_pairs) != len(lower_slopes) + 1:
        raise ValueError("Canonical input must have p upper and p-1 lower parameters")
    if not isinstance(eps, s.Symbol) or not isinstance(z, s.Expr) or z.has(eps):
        raise ValueError("Use a symbolic regulator and an epsilon-independent SymPy argument")
    upper = []
    for pair in upper_pairs:
        if not isinstance(pair, (tuple, list, s.Tuple)) or len(pair) != 2:
            raise ValueError("Each upper parameter must be a (base, slope) pair")
        base, slope = _rational(pair[0]), _slope(pair[1], eps)
        if base not in (0, 1):
            raise ValueError("Canonical upper bases must be 0 or 1")
        upper.append((base, slope))
    lower = [_slope(value, eps) for value in lower_slopes]
    zeros = [slope for base, slope in upper if base == 0]
    count = len(zeros)
    amplitude = s.prod(zeros)
    if amplitude == 0 or count > order:
        return s.S.One
    residual_order = order - count
    harmonic_log = {}
    rational_log = {}
    for weight in range(1, residual_order + 1):
        prefactor = s.Rational((-1) ** (weight + 1), weight)
        lower_power = sum((slope**weight for slope in lower), s.S.Zero)
        harmonic_log[weight] = {
            (weight,): prefactor
            * (sum((slope**weight for _, slope in upper), s.S.Zero) - lower_power)
        }
        rational_log[weight] = prefactor * (
            sum((slope**weight for base, slope in upper if base == 1), s.S.Zero) - lower_power
        )
    harmonic_coefficients = harmonic_exp(harmonic_log, residual_order)
    rational_coefficients = [s.S.One]
    for degree in range(1, residual_order + 1):
        rational_coefficients.append(
            sum(
                (
                    weight * rational_log[weight] * rational_coefficients[degree - weight]
                    for weight in range(1, degree + 1)
                ),
                s.S.Zero,
            )
            / degree
        )
    result = s.S.One
    for degree in range(residual_order + 1):
        for rational_degree in range(degree + 1):
            for word, coefficient in harmonic_coefficients[degree - rational_degree].items():
                result += (
                    amplitude
                    * eps ** (count + degree)
                    * coefficient
                    * rational_coefficients[rational_degree]
                    * generating_sum(count + rational_degree, word, z)
                )
    return s.expand(result)
