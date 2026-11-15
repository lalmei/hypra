"""Restrained SymPy nseries composition for exact HPL zero germs."""

import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.symbolic import HPL

Z = s.Symbol("z")


def test_composite_sqrt_argument_has_exact_puiseux_coefficients():
    actual = s.series(HPL((0, 1), s.sqrt(Z)), Z, 0, 3)
    expected = (
        s.sqrt(Z) + Z / 4 + Z ** s.Rational(3, 2) / 9 + Z**2 / 16 + Z ** s.Rational(5, 2) / 25
    )
    assert actual.removeO() == expected
    assert actual.getO() == s.Order(Z**3, Z)


def test_rational_composition_agrees_with_elementary_identity():
    argument = Z / (1 + Z)
    actual = s.series(HPL((1,), argument), Z, 0, 5)
    expected = s.series(-s.log(1 - argument), Z, 0, 5)
    assert s.expand(actual - expected) == s.Order(Z**5, Z)


def test_negative_argument_analytic_hpl_and_unknown_sign_rejection():
    actual = s.series(HPL((-1,), -s.sqrt(Z)), Z, 0, 3)
    expected = s.series(s.log(1 - s.sqrt(Z)), Z, 0, 3)
    assert s.expand(actual - expected) == s.Order(Z**3, Z)
    a = s.Symbol("a")
    with pytest.raises(UnsupportedExpressionError, match="known nonzero"):
        s.series(HPL((1,), a * Z), Z, 0, 3)


def test_leading_term_avoids_discrete_word_derivative():
    assert HPL((0, 1), s.sqrt(Z)).as_leading_term(Z) == s.sqrt(Z)
    assert HPL((1, 1), Z).as_leading_term(Z) == Z**2 / 2
    assert HPL((0, -1, 1), Z**2).as_leading_term(Z) == Z**4 / 4


def test_empty_and_all_zero_words_are_exact():
    assert s.series(HPL((), s.sqrt(Z)), Z, 0, 3) == 1
    actual = s.series(HPL((0, 0), s.sqrt(Z)), Z, 0, 3)
    assert actual.removeO() == s.log(Z) ** 2 / 8
    assert HPL((0,), Z).as_leading_term(Z) == s.log(Z)


def test_constant_hpl_is_an_unchanged_series_coefficient():
    constant = HPL((0, 1), s.Rational(1, 2))
    assert s.series(constant, Z, 0, 3) == constant


@pytest.mark.parametrize("argument", [1 + Z, 1 / Z, s.exp(-1 / Z)])
def test_nonzero_or_unknown_germ_fails_clearly(argument):
    with pytest.raises(UnsupportedExpressionError):
        s.series(HPL((0, 1), argument), Z, 0, 3)


def test_mixed_log_germ_directed_to_explicit_series_api():
    with pytest.raises(UnsupportedExpressionError, match="hpl_series"):
        s.series(HPL((1, 0), s.sqrt(Z)), Z, 0, 3)
