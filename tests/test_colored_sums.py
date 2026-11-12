"""Exact defining-sum and Pochhammer checks, independent of the HPL color map."""

from functools import cache

import pytest
import sympy as s

from hypra.colored_sums import (
    canonical_colored_hyper_series,
    colored_add,
    colored_exp,
    colored_generating_sum,
    colored_product,
    colored_quasi_shuffle,
    colored_scale,
    colored_word_to_hpl,
    even_generating_sum,
)
from hypra.symbolic import HPL

E, Z = s.symbols("eps z")
T = s.Symbol("t", positive=True)
HALF = s.Rational(1, 2)


@cache
def finite_colored_sum(n, word):
    if not word:
        return s.S.One
    index, sign = word[0]
    return sum(
        (
            sign**k * finite_colored_sum(k - 1, word[1:]) / s.Integer(k) ** index
            for k in range(1, n + 1)
        ),
        s.S.Zero,
    )


def evaluate_harmonic(poly, n):
    return sum((value * finite_colored_sum(n, word) for word, value in poly.items()), s.S.Zero)


@cache
def hpl_taylor_coefficients(word, nmax):
    """Independent integration of the HPL kernels at zero, not a color inversion."""
    if not word:
        return (s.S.One,) + (s.S.Zero,) * nmax
    assert word[-1] != 0, "Only HPLs analytic at zero are used in this test"
    tail = hpl_taylor_coefficients(word[1:], nmax)
    coefficients = [s.S.Zero]
    for n in range(1, nmax + 1):
        if word[0] == 0:
            value = tail[n] / n
        elif word[0] == 1:
            value = sum(tail[:n], s.S.Zero) / n
        elif word[0] == -1:
            value = sum(((-1) ** (n - 1 - k) * tail[k] for k in range(n)), s.S.Zero) / n
        else:
            raise AssertionError("Not an HPL letter")
        coefficients.append(value)
    return tuple(coefficients)


def replace_hpl_by_taylor(expression, nmax):
    substitutions = {}
    for atom in expression.atoms(HPL):
        coefficients = hpl_taylor_coefficients(tuple(atom.word), nmax)
        substitutions[atom] = sum(value * T**k for k, value in enumerate(coefficients))
    return expression.xreplace(substitutions)


def test_colored_merge_signs_and_multiplicities():
    word = ((1, -1),)
    expected = {((1, -1), (1, -1)): 2, ((2, 1),): 1}
    assert colored_quasi_shuffle(word, word) == expected
    assert colored_quasi_shuffle((), word) == {word: 1}
    result = colored_quasi_shuffle(word, word)
    result.clear()
    assert colored_quasi_shuffle(word, word) == expected


@pytest.mark.parametrize(
    "left,right",
    [
        (((1, -1),), ((1, 1),)),
        (((2, -1), (1, 1)), ((1, -1), (2, -1))),
        (((1, 1),), ((1, -1), (1, -1))),
    ],
)
def test_colored_quasi_shuffle_matches_literal_sums(left, right):
    product = colored_quasi_shuffle(left, right)
    for n in range(1, 6):
        assert evaluate_harmonic(product, n) == finite_colored_sum(n, left) * finite_colored_sum(
            n, right
        )


def test_colored_algebra_and_exponential():
    one = {(): s.S.One}
    poly = {((1, -1),): s.Integer(2), ((2, 1),): s.Rational(1, 3)}
    assert colored_product(poly, one) == poly
    assert colored_add(poly, colored_scale(poly, -1)) == {}
    coefficients = colored_exp({1: poly}, 3)
    for n in range(1, 5):
        for degree, coefficient in enumerate(coefficients):
            assert evaluate_harmonic(coefficient, n) == evaluate_harmonic(
                poly, n
            ) ** degree / s.factorial(degree)


def test_color_to_hpl_sign_conventions():
    assert colored_word_to_hpl(((1, -1),), T) == -HPL((-1,), T)
    assert colored_word_to_hpl(((2, -1), (1, -1)), T) == -HPL((0, -1, 1), T)
    assert colored_word_to_hpl(((1, 1), (2, -1)), T) == -HPL((1, 0, -1), T)


@pytest.mark.parametrize(
    "power,word",
    [
        (0, ()),
        (0, ((1, -1),)),
        (0, ((2, -1), (1, 1))),
        (1, ()),
        (2, ((1, -1),)),
        (1, ((2, -1), (1, -1))),
    ],
)
@pytest.mark.parametrize("delta", [-1, 1])
def test_generating_maps_by_independent_hpl_kernels(power, word, delta):
    expr = colored_generating_sum(power, word, T, delta)
    series = s.series(replace_hpl_by_taylor(expr, 7), T, 0, 8).removeO().expand()
    for n in range(1, 8):
        expected = delta**n * finite_colored_sum(n - 1, word) / s.Integer(n) ** power
        assert series.coeff(T, n) == expected


def test_even_projection_removes_odd_coefficients():
    word = ((1, -1), (1, 1))
    expr = even_generating_sum(1, word, T)
    series = replace_hpl_by_taylor(expr, 8).expand()
    for n in range(1, 9):
        expected = finite_colored_sum(n - 1, word) / n if n % 2 == 0 else 0
        assert series.coeff(T, n) == expected


@pytest.mark.parametrize(
    "upper,lower,order",
    [
        ([(1, 2), (HALF, -1)], [(HALF, 1)], 3),
        ([(0, -1), (HALF, 2)], [(HALF, -2)], 4),
        ([(1, 0), (0, 1), (HALF, s.Rational(1, 3))], [(1, 2), (HALF, -1)], 3),
        ([(0, 1), (0, -2), (HALF, 1)], [(1, 2), (HALF, 3)], 3),
        ([(1, 1), (HALF, 2), (HALF, -1)], [(HALF, 1), (HALF, -2)], 3),
    ],
)
def test_matched_half_series_against_direct_pochhammer_coefficients(upper, lower, order):
    nmax = 3
    computed = canonical_colored_hyper_series(upper, lower, Z, E, order)
    computed = replace_hpl_by_taylor(computed, 2 * nmax).subs(Z, T**2)
    computed = s.series(computed, T, 0, 2 * nmax + 1).removeO().expand()
    assert s.simplify(computed.coeff(T, 0) - 1) == 0
    for n in range(1, nmax + 1):
        assert s.expand(computed.coeff(T, 2 * n - 1)) == 0
        numerator = s.Poly(s.prod(s.rf(base + slope * E, n) for base, slope in upper), E)
        denominator = s.Poly(
            s.factorial(n) * s.prod(s.rf(base + slope * E, n) for base, slope in lower), E
        )
        coefficients = []
        for degree in range(order + 1):
            coefficients.append(
                (
                    numerator.nth(degree)
                    - sum(
                        denominator.nth(k) * coefficients[degree - k] for k in range(1, degree + 1)
                    )
                )
                / denominator.nth(0)
            )
        expected = sum(value * E**degree for degree, value in enumerate(coefficients))
        assert s.expand(computed.coeff(T, 2 * n) - expected) == 0


def test_cancellation_zero_upper_and_integer_fallback():
    assert canonical_colored_hyper_series([(0, 0), (HALF, 1)], [(HALF, 2)], Z, E, 3) == 1
    assert (
        canonical_colored_hyper_series([(0, 1), (0, 2), (HALF, 1)], [(1, 3), (HALF, 2)], Z, E, 1)
        == 1
    )
    integer = canonical_colored_hyper_series([(0, -1), (1, 0)], [(1, -2)], Z, E, 1)
    assert integer == 1 - E * HPL((1,), Z)


@pytest.mark.parametrize(
    "call",
    [
        lambda: colored_quasi_shuffle(((0, -1),), ()),
        lambda: colored_quasi_shuffle(((1, 0),), ()),
        lambda: colored_quasi_shuffle(((1, True),), ()),
        lambda: colored_quasi_shuffle(((1.0, 1),), ()),
        lambda: colored_generating_sum(1, (), T, 0),
        lambda: canonical_colored_hyper_series([(1, 0), (HALF, 1)], [(1, 1)], Z, E, 2),
        lambda: canonical_colored_hyper_series([(1, 0), (HALF, 1.0)], [(HALF, 2)], Z, E, 2),
        lambda: canonical_colored_hyper_series(
            [(1, 0), (HALF, 1)], [(HALF, 2)], s.Float(".5"), E, 2
        ),
    ],
)
def test_invalid_colored_inputs_fail_closed(call):
    with pytest.raises((ValueError, TypeError)):
        call()


def test_internal_colored_formal_slopes_match_pochhammer():
    delta = s.Symbol("delta")
    upper, lower = [(0, delta), (HALF, 1 + delta)], [(HALF, 2 - delta)]
    computed = canonical_colored_hyper_series(upper, lower, Z, E, 2)
    computed = replace_hpl_by_taylor(computed, 4).subs(Z, T**2)
    computed = s.series(computed, T, 0, 5).removeO().expand()
    for n in (1, 2):
        numerator = s.Poly(s.prod(s.rf(base + slope * E, n) for base, slope in upper), E)
        denominator = s.Poly(
            s.factorial(n) * s.prod(s.rf(base + slope * E, n) for base, slope in lower), E
        )
        coefficients = []
        for degree in range(3):
            coefficients.append(
                (
                    numerator.nth(degree)
                    - sum(
                        denominator.nth(k) * coefficients[degree - k] for k in range(1, degree + 1)
                    )
                )
                / denominator.nth(0)
            )
        expected = sum(value * E**degree for degree, value in enumerate(coefficients))
        assert s.expand(computed.coeff(T, 2 * n) - expected) == 0


def test_colored_regulator_dependent_or_nonpolynomial_slopes_rejected():
    delta = s.Symbol("delta")
    for slope in (E + delta, 1 / delta):
        with pytest.raises((ValueError, TypeError)):
            canonical_colored_hyper_series([(0, slope), (HALF, 1)], [(HALF, 2)], Z, E, 2)
