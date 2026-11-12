"""Matched-half canonical hypergeometric expansions by colored strict sums.

Color convention: H_N((a1,s1),...) sums prod(sj**nj/nj**aj) over
N >= n1 > ... > nd >= 1, with signs sj in {-1,1}. In particular the
negative color is (-1)**n, not (-1)**(n-1). Products merge indices by
addition and colors by multiplication, with a positive stuffle sign.

For equally many half-integer upper and lower bases, their unperturbed
Pochhammer factors cancel. Doubling the sum index and projecting onto even
indices then reduces all coefficient sums to the HPL alphabet {0,1,-1}.
This implements the defining-product derivation, not imported expansion tables.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from functools import lru_cache

import sympy as s

from .nested_sums import _coefficient, _order, _rational, _slope, canonical_hyper_series
from .symbolic import HPL

ColoredWord = tuple[tuple[int, int], ...]
ColoredPolynomial = dict[ColoredWord, s.Expr]


def _word(value) -> ColoredWord:
    if not isinstance(value, (tuple, list, s.Tuple)):
        raise TypeError("A colored word must be a sequence of (positive index, sign) pairs")
    result = []
    for pair in value:
        if not isinstance(pair, (tuple, list, s.Tuple)) or len(pair) != 2:
            raise ValueError("A colored letter must be an (index, sign) pair")
        index, sign = pair
        if isinstance(index, bool) or not isinstance(index, (int, s.Integer)) or index <= 0:
            raise ValueError("Colored indices must be positive integers")
        if isinstance(sign, bool) or not isinstance(sign, (int, s.Integer)) or sign not in (-1, 1):
            raise ValueError("Colored signs must be exact -1 or +1")
        result.append((int(index), int(sign)))
    return tuple(result)


def _clean(poly: Mapping) -> ColoredPolynomial:
    if not isinstance(poly, Mapping):
        raise TypeError("Colored polynomials must map words to exact coefficients")
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
def _quasi_shuffle(left: ColoredWord, right: ColoredWord) -> tuple[tuple[ColoredWord, int], ...]:
    if not left:
        return ((right, 1),)
    if not right:
        return ((left, 1),)
    terms = {}
    merged = (left[0][0] + right[0][0], left[0][1] * right[0][1])
    for prefix, rest_left, rest_right in (
        (left[0], left[1:], right),
        (right[0], left, right[1:]),
        (merged, left[1:], right[1:]),
    ):
        for tail, multiplicity in _quasi_shuffle(rest_left, rest_right):
            word = (prefix, *tail)
            terms[word] = terms.get(word, 0) + multiplicity
    return tuple(sorted(terms.items()))


def colored_quasi_shuffle(left, right) -> dict[ColoredWord, int]:
    """Strict colored stuffle product, including coincident-word multiplicities."""
    return dict(_quasi_shuffle(_word(left), _word(right)))


def colored_add(*polynomials: Mapping) -> ColoredPolynomial:
    result = {}
    for poly in polynomials:
        for word, coefficient in _clean(poly).items():
            result[word] = result.get(word, s.S.Zero) + coefficient
    return _clean(result)


def colored_scale(polynomial: Mapping, coefficient) -> ColoredPolynomial:
    scalar = _coefficient(coefficient)
    return _clean({word: scalar * value for word, value in _clean(polynomial).items()})


def colored_product(left: Mapping, right: Mapping) -> ColoredPolynomial:
    terms = {}
    for left_word, left_coefficient in _clean(left).items():
        for right_word, right_coefficient in _clean(right).items():
            for word, multiplicity in _quasi_shuffle(left_word, right_word):
                terms[word] = terms.get(word, s.S.Zero) + (
                    multiplicity * left_coefficient * right_coefficient
                )
    return _clean(terms)


def colored_exp(log_coefficients: Mapping[int, Mapping], order: int) -> list[ColoredPolynomial]:
    """Formal exp coefficients via n F_n = sum_m m L_m * F_(n-m)."""
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
        for weight in range(1, degree + 1):
            total = colored_add(
                total,
                colored_scale(
                    colored_product(logarithm.get(weight, {}), coefficients[degree - weight]),
                    weight,
                ),
            )
        coefficients.append(colored_scale(total, s.Rational(1, degree)))
    return coefficients


def colored_word_to_hpl(word, argument):
    """Map a colored multiple polylogarithm to expanded HPL letters.

    The input denotes sum_{n1>...>nd>0} argument**n1 * prod(sign_j**nj/nj**aj).
    Nonzero HPL letters are cumulative products of these signs. Multiplication
    by the product of those cumulative signs converts the standard HPL kernels
    to the chosen colored-sum convention. The empty word is one.
    """
    word = _word(word)
    if not word:
        return s.S.One
    letters = []
    cumulative, prefactor = 1, 1
    for index, sign in word:
        cumulative *= sign
        prefactor *= cumulative
        letters.extend([0] * (index - 1))
        letters.append(cumulative)
    return prefactor * HPL(tuple(letters), argument)


def colored_generating_sum(power: int, word, argument, outer_color: int = 1):
    """Sum delta**n*x**n/n**power H_(n-1)(word), initially |x|<1.

    ``delta`` is outer_color. For power zero the geometric sum must be done
    before HPL conversion; simply adding an index zero would be incorrect.
    """
    if type(power) is not int or power < 0:
        raise ValueError("Generating-sum power must be a nonnegative integer")
    if type(outer_color) is not int or outer_color not in (-1, 1):
        raise ValueError("Outer color must be -1 or +1")
    word = _word(word)
    if power:
        return colored_word_to_hpl(((power, outer_color), *word), argument)
    prefactor = outer_color * argument / (1 - outer_color * argument)
    if not word:
        return prefactor
    first, *tail = word
    return prefactor * colored_word_to_hpl(((first[0], outer_color * first[1]), *tail), argument)


def even_generating_sum(power: int, word, argument):
    """Project the generating series to even outer indices, at argument t."""
    return (
        colored_generating_sum(power, word, argument, 1)
        + colored_generating_sum(power, word, argument, -1)
    ) / 2


def _pairs(pairs, allowed, label, regulator=None):
    if not isinstance(pairs, Sequence):
        raise TypeError(f"{label} parameters must be a finite sequence")
    result = []
    for pair in pairs:
        if not isinstance(pair, (tuple, list, s.Tuple)) or len(pair) != 2:
            raise ValueError(f"Each {label} parameter must be a (base, slope) pair")
        base, slope = _rational(pair[0]), _slope(pair[1], regulator)
        if base not in allowed:
            raise ValueError(f"Unsupported canonical {label} parameter base")
        result.append((base, slope))
    return result


def canonical_colored_hyper_series(upper_pairs, lower_pairs, z, eps, order):
    """Exact canonical matched-half pF(p-1) series as unreduced HPLs at sqrt(z).

    Upper bases may be 0, 1 or 1/2; lower bases 1 or 1/2. Half-base counts must
    match. Internal slopes may be exact algebraic polynomials in formal
    regularization symbols; public planning still restricts slopes to rationals.
    The caller enforces the real domain
    0<z<1 before using this expansion, and owns the truncation/prefactor budget.
    This routine does not attempt branch transformations or emit an O term.

    With integer-base slope power differences U_m (all integer uppers minus
    integer lowers), V_m (base-one uppers minus integer lowers), and half-base
    difference D_m, the doubled-index logarithm has harmonic coefficients
    (-1)**(m+1)*2**(m-1)/m times [(U_m+D_m)H^+ + (U_m-D_m)H^-],
    plus (-1)**(m+1)*2**m*V_m/(m*(2n)**m). The global coefficient is
    A*eps**r*2**r/(2n)**r, where r counts zero bases and A is their slope product.
    """
    _order(order)
    half = s.Rational(1, 2)
    upper = _pairs(upper_pairs, (0, 1, half), "upper", eps)
    lower = _pairs(lower_pairs, (1, half), "lower", eps)
    if len(upper) != len(lower) + 1:
        raise ValueError("Canonical colored input must be pF(p-1)")
    if not isinstance(eps, s.Symbol) or not isinstance(z, s.Expr) or z.has(eps):
        raise ValueError("Use a symbolic regulator and an epsilon-independent SymPy argument")
    if z.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise ValueError("Canonical arguments must be exact and finite")
    upper_half = [slope for base, slope in upper if base == half]
    lower_half = [slope for base, slope in lower if base == half]
    if len(upper_half) != len(lower_half):
        raise ValueError("Canonical colored expansion requires matched half-base counts")
    if not upper_half:
        return canonical_hyper_series(upper, [slope for _, slope in lower], z, eps, order)
    zeros = [slope for base, slope in upper if base == 0]
    count, amplitude = len(zeros), s.prod(zeros)
    if amplitude == 0 or count > order:
        return s.S.One
    residual_order = order - count
    upper_integer = [slope for base, slope in upper if base != half]
    upper_one = [slope for base, slope in upper if base == 1]
    lower_integer = [slope for base, slope in lower if base == 1]
    harmonic_log, rational_log = {}, {}
    for weight in range(1, residual_order + 1):
        prefactor = s.Rational((-1) ** (weight + 1), weight)
        lower_power = sum((slope**weight for slope in lower_integer), s.S.Zero)
        integer_difference = sum((slope**weight for slope in upper_integer), s.S.Zero) - lower_power
        one_difference = sum((slope**weight for slope in upper_one), s.S.Zero) - lower_power
        half_difference = sum((slope**weight for slope in upper_half), s.S.Zero) - sum(
            (slope**weight for slope in lower_half), s.S.Zero
        )
        harmonic_log[weight] = {
            ((weight, 1),): prefactor * 2 ** (weight - 1) * (integer_difference + half_difference),
            ((weight, -1),): prefactor * 2 ** (weight - 1) * (integer_difference - half_difference),
        }
        rational_log[weight] = prefactor * 2**weight * one_difference
    harmonic_coefficients = colored_exp(harmonic_log, residual_order)
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
    t = s.sqrt(z)
    result = s.S.One
    for degree in range(residual_order + 1):
        for rational_degree in range(degree + 1):
            for word, coefficient in harmonic_coefficients[degree - rational_degree].items():
                result += (
                    amplitude
                    * 2**count
                    * eps ** (count + degree)
                    * coefficient
                    * rational_coefficients[rational_degree]
                    * even_generating_sum(count + rational_degree, word, t)
                )
    return s.expand(result)
