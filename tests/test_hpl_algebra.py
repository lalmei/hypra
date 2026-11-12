"""Exact integration/shuffle tests independent of hypergeometric recurrences."""

import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import (
    integrate_from_zero,
    inverse_theta_plus,
    shuffle_expand,
    shuffle_product,
    shuffle_words,
)
from hypra.symbolic import HPL, reduce_known_hpl

Z = s.Symbol("z", positive=True)


def h(*word):
    return HPL(s.Tuple(*word), Z)


def assert_derivative(primitive, integrand):
    residual = shuffle_expand(s.diff(primitive, Z) - integrand, Z)
    assert s.cancel(residual) == 0


def test_shuffle_multiplicity_and_empty_word():
    assert shuffle_words((1,), (1,)) == {(1, 1): 2}
    assert shuffle_words((), (0, 1)) == {(0, 1): 1}
    assert shuffle_words((0,), (1,)) == {(0, 1): 1, (1, 0): 1}
    assert shuffle_product(h(1), h(1)) == 2 * h(1, 1)
    assert shuffle_expand(h(1) ** 3, Z) == 6 * h(1, 1, 1)


def test_log_and_polylog_conversion_preserves_signs():
    expression = s.log(Z) * s.log(1 - Z) + s.polylog(3, -Z)
    assert s.expand(shuffle_expand(expression, Z) + h(0, 1) + h(1, 0) + h(0, 0, -1)) == 0
    assert shuffle_expand(s.log(1 + Z), Z) == h(-1)


def test_polynomial_and_simple_kernels():
    assert integrate_from_zero(3 * Z**2 + 2, Z) == Z**3 + 2 * Z
    assert integrate_from_zero(h(1) / Z, Z) == h(0, 1)
    assert integrate_from_zero(h(0, 1) / (1 + Z), Z) == h(-1, 0, 1)
    assert integrate_from_zero(1 / (1 - Z), Z) == h(1)


@pytest.mark.parametrize(
    "expression",
    [
        Z**3 * h(0, 1),
        h(1) / (1 - Z) ** 3,
        (Z + 2) * h(-1, 0, 1) / (1 + Z) ** 2,
        h(0) * h(1) * (1 + Z),
        s.log(Z) ** 2 * h(0, 1) / Z,
    ],
)
def test_integrated_primitive_differentiates_exactly(expression):
    primitive = integrate_from_zero(expression, Z)
    assert not primitive.has(s.Integral)
    assert_derivative(primitive, expression)


def test_logarithmic_origin_integrable_case():
    primitive = integrate_from_zero(s.log(Z), Z)
    assert s.expand(reduce_known_hpl(primitive) - (Z * s.log(Z) - Z)) == 0


def test_cancelled_divergences_are_combined_before_origin_test():
    expression = (h(1) - Z) / Z**2
    primitive = integrate_from_zero(expression, Z)
    assert_derivative(primitive, expression)
    assert s.limit(reduce_known_hpl(primitive), Z, 0, dir="+") == 0
    assert s.expand(primitive - (1 + h(1) - h(1) / Z)) == 0


def test_repeated_origin_pole_finite_cancellation():
    expression = (h(1) - Z - Z**2 / 2) / Z**3
    primitive = integrate_from_zero(expression, Z)
    assert_derivative(primitive, expression)
    assert s.limit(reduce_known_hpl(primitive), Z, 0, dir="+") == 0


def test_parameter_coefficients_are_exact():
    a = s.Symbol("a")
    primitive = integrate_from_zero(a * h(1) / Z + s.zeta(3), Z)
    assert primitive == a * h(0, 1) + Z * s.zeta(3)


@pytest.mark.parametrize("expression", [1 / Z, h(1) / Z**2, s.log(Z) / Z, Z**-2])
def test_divergent_origin_is_rejected(expression):
    with pytest.raises(UnsupportedExpressionError, match="origin"):
        integrate_from_zero(expression, Z)


@pytest.mark.parametrize("expression", [1 / (2 - Z), s.sin(Z), 1 / h(1), s.log(2 - Z)])
def test_unsupported_inputs_fail_closed(expression):
    with pytest.raises(UnsupportedExpressionError):
        integrate_from_zero(expression, Z)


def test_float_and_wrong_argument_rejected():
    with pytest.raises(UnsupportedExpressionError):
        integrate_from_zero(s.Float("0.1") * Z, Z)
    with pytest.raises(UnsupportedExpressionError):
        integrate_from_zero(HPL(s.Tuple(1), 1 - Z), Z)


@pytest.mark.parametrize("n", [1, 2, 4])
def test_inverse_theta_selects_regular_solution(n):
    integrand = h(0, 1) + Z * h(1)
    solution = inverse_theta_plus(integrand, Z, n)
    assert s.cancel(shuffle_expand(Z * s.diff(solution, Z) + n * solution - integrand, Z)) == 0
    assert s.limit(reduce_known_hpl(solution), Z, 0, dir="+") == 0


def test_resource_limits_fail_closed():
    with pytest.raises(UnsupportedExpressionError):
        shuffle_words((1,) * 25, ())
    with pytest.raises(UnsupportedExpressionError):
        integrate_from_zero(1 / (1 - Z) ** 33, Z)
    with pytest.raises(UnsupportedExpressionError):
        inverse_theta_plus(Z, Z, 0)


def test_hpl_local_series_matches_log_and_dilog_series():
    from hypra.hpl_algebra import hpl_series

    expression = h(0, 1) + h(1, 1) + h(-1)
    actual = hpl_series(expression, Z, 6)
    expected = s.series(s.polylog(2, Z) + s.log(1 - Z) ** 2 / 2 + s.log(1 + Z), Z, 0, 6).removeO()
    assert s.expand(actual - expected) == 0


def test_hpl_local_series_handles_trailing_zero_and_rational_pole():
    from hypra.hpl_algebra import hpl_series

    # H_(1,0)= -log(z)log(1-z)-Li2(z), with vanishing regularized origin.
    expression = h(1, 0) / Z**2
    actual = hpl_series(expression, Z, 4)
    expected = s.series((-s.log(Z) * s.log(1 - Z) - s.polylog(2, Z)) / Z**2, Z, 0, 4).removeO()
    assert s.expand(actual - expected) == 0


def test_pre_expansion_budget_rejects_huge_multinomial():
    with pytest.raises(UnsupportedExpressionError, match="budget"):
        shuffle_expand((1 + h(1) + h(0)) ** 100, Z)
