"""Independent finite-sum tests for canonical hypergeometric coefficient algebra."""

from functools import cache

import pytest
import sympy as s

from hypra.nested_sums import (
    canonical_hyper_series,
    generating_sum,
    harmonic_add,
    harmonic_exp,
    harmonic_product,
    harmonic_scale,
    harmonic_word_to_hpl,
    quasi_shuffle,
)
from hypra.symbolic import HPL

E, Z = s.symbols("eps z")


@cache
def finite_sum(n, word):
    """Strict nested sum by its definition, independent of stuffle implementation."""
    if not word:
        return s.S.One
    return sum(
        (finite_sum(k - 1, word[1:]) / s.Integer(k) ** word[0] for k in range(1, n + 1)), s.S.Zero
    )


def evaluate_harmonic(poly, n):
    return sum((value * finite_sum(n, word) for word, value in poly.items()), s.S.Zero)


def indices_from_hpl(atom):
    indices, zeros = [], 0
    for letter in atom.word:
        if letter == 0:
            zeros += 1
        elif letter == 1:
            indices.append(zeros + 1)
            zeros = 0
        else:
            raise AssertionError("Canonical sums use only letters 0 and 1")
    assert zeros == 0
    return tuple(indices)


def hpl_taylor(atom, nmax):
    word = indices_from_hpl(atom)
    if not word:
        return s.S.One
    return sum(
        (Z**n * finite_sum(n - 1, word[1:]) / s.Integer(n) ** word[0] for n in range(1, nmax + 1)),
        s.S.Zero,
    )


def test_stuffle_multiplicity_and_merge_sign():
    assert quasi_shuffle((1,), (1,)) == {(1, 1): 2, (2,): 1}
    assert quasi_shuffle((1,), (1, 1)) == {(1, 1, 1): 3, (1, 2): 1, (2, 1): 1}
    assert quasi_shuffle((), (2, 1)) == {(2, 1): 1}
    value = quasi_shuffle((1,), (1,))
    value[(2,)] = -100
    assert quasi_shuffle((1,), (1,))[(2,)] == 1


@pytest.mark.parametrize("left,right", [((1,), (2,)), ((1, 2), (2, 1)), ((2,), (1, 1)), ((), (3,))])
def test_quasi_shuffle_matches_finite_nested_sums(left, right):
    product = quasi_shuffle(left, right)
    for n in range(1, 7):
        assert finite_sum(n, left) * finite_sum(n, right) == evaluate_harmonic(product, n)


def test_algebra_and_formal_exponential():
    one = {(): s.S.One}
    h1 = {(1,): s.S.One}
    assert harmonic_product(one, h1) == h1
    assert harmonic_add(h1, harmonic_scale(h1, -1)) == {}
    coefficients = harmonic_exp({1: h1}, 3)
    assert coefficients == [
        one,
        h1,
        {(1, 1): s.S.One, (2,): s.Rational(1, 2)},
        {
            (1, 1, 1): s.S.One,
            (1, 2): s.Rational(1, 2),
            (2, 1): s.Rational(1, 2),
            (3,): s.Rational(1, 6),
        },
    ]
    for n in range(1, 5):
        for degree, coefficient in enumerate(coefficients):
            assert evaluate_harmonic(coefficient, n) == finite_sum(n, (1,)) ** degree / s.factorial(
                degree
            )


def test_hpl_and_geometric_sum_maps():
    assert harmonic_word_to_hpl((2, 1, 3), Z) == HPL((0, 1, 1, 0, 0, 1), Z)
    assert generating_sum(2, (1,), Z) == HPL((0, 1, 1), Z)
    assert generating_sum(0, (), Z) == Z / (1 - Z)
    assert generating_sum(0, (2,), Z) == Z / (1 - Z) * HPL((0, 1), Z)
    for power, word in [(0, ()), (0, (1, 2)), (1, (1,)), (3, (2,))]:
        expr = generating_sum(power, word, Z)
        expr = expr.xreplace({h: hpl_taylor(h, 6) for h in expr.atoms(HPL)})
        expr = s.series(expr, Z, 0, 7).removeO().expand()
        for n in range(1, 7):
            assert expr.coeff(Z, n) == finite_sum(n - 1, word) / s.Integer(n) ** power


@pytest.mark.parametrize(
    "upper,lower",
    [
        ([(1, 2)], []),
        ([(0, s.Rational(-1, 2))], []),
        ([(1, 0), (0, -1)], [-2]),
        ([(1, 2), (0, s.Rational(1, 2))], [3]),
        ([(0, 1), (0, -2), (1, s.Rational(1, 3))], [2, -1]),
        ([(1, 1), (1, -1), (1, 2)], [s.Rational(1, 2), -2]),
        ([(0, 1), (0, -1), (0, 2), (1, 3)], [1, 2, -1]),
    ],
)
def test_canonical_series_against_direct_pochhammer_coefficients(upper, lower):
    order, nmax = 4, 4
    computed = canonical_hyper_series(upper, lower, Z, E, order)
    computed = computed.xreplace({h: hpl_taylor(h, nmax) for h in computed.atoms(HPL)})
    computed = s.series(computed, Z, 0, nmax + 1).removeO().expand()
    assert s.simplify(computed.coeff(Z, 0) - 1) == 0
    for n in range(1, nmax + 1):
        numerator = s.Poly(s.prod(s.rf(base + slope * E, n) for base, slope in upper), E)
        denominator = s.Poly(s.factorial(n) * s.prod(s.rf(1 + slope * E, n) for slope in lower), E)
        # Direct finite-product rational division, independent of logarithmic
        # harmonic expansion. Avoid expensive general-purpose series inversion.
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
        assert s.expand(computed.coeff(Z, n) - expected) == 0


def test_zero_upper_and_late_start_are_exact():
    assert canonical_hyper_series([(0, 0), (1, 1)], [1], Z, E, 4) == 1
    assert canonical_hyper_series([(0, 1), (0, 2)], [1], Z, E, 1) == 1
    assert s.simplify(canonical_hyper_series([(1, 0)], [], Z, E, 0) - 1 / (1 - Z)) == 0


@pytest.mark.parametrize(
    "call",
    [
        lambda: quasi_shuffle((0,), (1,)),
        lambda: quasi_shuffle((-1,), (1,)),
        lambda: quasi_shuffle((1.0,), (1,)),
        lambda: harmonic_scale({(1,): s.S.One}, 0.1),
        lambda: harmonic_exp({0: {(): s.S.One}}, 1),
        lambda: generating_sum(-1, (), Z),
        lambda: canonical_hyper_series([(2, 1)], [], Z, E, 2),
        lambda: canonical_hyper_series([(1, 0.1)], [], Z, E, 2),
        lambda: canonical_hyper_series([(1, 1)], [1], Z, E, 2),
        lambda: canonical_hyper_series([(1, 1)], [], Z + E, E, 2),
    ],
)
def test_malformed_inputs_rejected(call):
    with pytest.raises((ValueError, TypeError)):
        call()


def test_internal_formal_polynomial_slopes_are_exact():
    delta = s.Symbol("delta")
    computed = canonical_hyper_series([(1, 0), (0, delta)], [2 - delta], Z, E, 2)
    expected = (
        1
        + delta * E * HPL((1,), Z)
        + E**2 * (delta * (delta - 2) * HPL((0, 1), Z) + 2 * delta * (delta - 1) * HPL((1, 1), Z))
    )
    assert s.expand(computed - expected) == 0
    algebraic = canonical_hyper_series([(0, s.sqrt(2) * delta)], [], Z, E, 1)
    assert algebraic == 1 + s.sqrt(2) * delta * E * HPL((1,), Z)


@pytest.mark.parametrize(
    "slope",
    [
        1 / s.Symbol("delta"),
        s.sin(s.Symbol("delta")),
        s.pi * s.Symbol("delta"),
        E + s.Symbol("delta"),
        s.Float("0.1"),
        s.oo,
    ],
)
def test_invalid_internal_slopes_fail_closed(slope):
    with pytest.raises((ValueError, TypeError)):
        canonical_hyper_series([(0, slope)], [], Z, E, 2)
