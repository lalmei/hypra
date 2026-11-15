"""Exact mathematical checks of contiguous shifts, independent of native code."""

import pytest
from sympy import Integer, Rational, factorial, integrate, rf, simplify, symbols


def coefficient(upper, lower, n):
    result = 1 / factorial(n)
    for a in upper:
        result *= rf(a, n)
    for b in lower:
        result /= rf(b, n)
    return result


@pytest.mark.parametrize("n", range(7))
def test_upper_raise_coefficient_identity(n):
    a, b, c = symbols("a b c", nonzero=True)
    old = coefficient([a, b], [c], n)
    new = coefficient([a + 1, b], [c], n)
    assert simplify(new - (n + a) / a * old) == 0


@pytest.mark.parametrize("n", range(7))
def test_lower_decrease_coefficient_identity(n):
    a, b, c = symbols("a b c", nonzero=True)
    old = coefficient([a, b], [c], n)
    new = coefficient([a, b], [c - 1], n)
    assert simplify(new - (n + c - 1) / (c - 1) * old) == 0


@pytest.mark.parametrize("n", range(7))
def test_lower_increase_coefficient_identity(n):
    a, b, c = symbols("a b c", nonzero=True)
    old = coefficient([a, b], [c], n)
    new = coefficient([a, b], [c + 1], n)
    assert simplify(new - c / (n + c) * old) == 0


@pytest.mark.parametrize("c", [1, 2, 3, Rational(3, 2)])
@pytest.mark.parametrize("n", range(4))
def test_lower_increase_integral_on_a_monomial(c, n):
    t, z = symbols("t z", positive=True)
    integrated = c * z ** (-c) * integrate(t ** (c - 1 + n), (t, 0, z))
    assert simplify(integrated - c / (c + n) * z**n) == 0


def test_lower_decrease_requires_one_deeper_order_even_if_final_cancels():
    eps = symbols("eps")
    exact_base_coefficient = 1 / (1 + eps)
    exact_shifted = simplify((1 + eps) / eps * exact_base_coefficient)
    shallow = exact_base_coefficient.series(eps, 0, 4).removeO()
    deep = exact_base_coefficient.series(eps, 0, 5).removeO()
    shallow_shifted = ((1 + eps) / eps * shallow).expand()
    deep_shifted = ((1 + eps) / eps * deep).expand()
    assert exact_shifted == 1 / eps
    assert shallow_shifted.coeff(eps, 3) == -1
    assert deep_shifted.coeff(eps, 3) == 0


@pytest.mark.parametrize("n", range(3, 8))
def test_negative_upper_tail_decomposition(n):
    eps = symbols("eps")
    upper = [-2 + eps, 1 + 2 * eps]
    lower = [2 + 3 * eps]
    shift = 3
    prefix_factor = coefficient(upper, lower, shift)
    tail_coefficient = coefficient(
        [Integer(1), *(a + shift for a in upper)],
        [Integer(shift + 1), *(b + shift for b in lower)],
        n - shift,
    )
    assert simplify(coefficient(upper, lower, n) - prefix_factor * tail_coefficient) == 0


@pytest.mark.parametrize("n", range(3, 8))
def test_exact_unit_upper_preserves_rank_in_tail(n):
    eps = symbols("eps")
    upper, lower = [Integer(1), -2 + eps], [2 + 3 * eps]
    shift = 3
    prefix_factor = coefficient(upper, lower, shift)
    reduced_tail = coefficient([Integer(1), 1 + eps], [5 + 3 * eps], n - shift)
    assert simplify(coefficient(upper, lower, n) - prefix_factor * reduced_tail) == 0


def test_exact_negative_upper_terminates_with_admissible_lower():
    for n in range(4, 9):
        assert coefficient([-3, Rational(1, 2)], [2], n) == 0


@pytest.mark.parametrize("integer_part", [0, -1, -2, -3])
def test_single_pochhammer_zero_threshold(integer_part):
    eps = symbols("eps")
    for n in range(7):
        polynomial = rf(integer_part + 2 * eps, n).expand()
        if n >= 1 - integer_part:
            assert polynomial.subs(eps, 0) == 0
            assert polynomial.coeff(eps, 1) != 0
        else:
            assert polynomial.subs(eps, 0) != 0


@pytest.mark.parametrize("n", range(4, 8))
def test_tail_identity_with_moving_negative_lower_parameters(n):
    eps = symbols("eps")
    upper = [-1 + eps, 1 + 2 * eps, 2 + 3 * eps]
    lower = [-2 + 4 * eps, -3 + 5 * eps]
    shift = 4
    factor = coefficient(upper, lower, shift)
    tail = coefficient(
        [Integer(1), *(a + shift for a in upper)],
        [Integer(shift + 1), *(b + shift for b in lower)],
        n - shift,
    )
    assert simplify(coefficient(upper, lower, n) - factor * tail) == 0


def test_tail_prefactor_pole_requires_deeper_tail_coefficients():
    eps = symbols("eps")
    upper = [-1 + eps, 1 + 2 * eps, 2 + 3 * eps]
    lower = [-2 + 4 * eps, -3 + 5 * eps]
    shift = 4
    factor = coefficient(upper, lower, shift)
    tail = coefficient(
        [Integer(1), *(a + shift for a in upper)],
        [Integer(shift + 1), *(b + shift for b in lower)],
        1,
    )
    assert factor.as_leading_term(eps).as_powers_dict()[eps] == -1
    target = 2
    shallow = tail.series(eps, 0, target + 1).removeO()
    deep = tail.series(eps, 0, target + 2).removeO()
    exact = coefficient(upper, lower, shift + 1)
    # Compare the omitted tail remainder before expansion of the complete factor.
    shallow_error = simplify(factor * (tail - shallow)).as_leading_term(eps)
    deep_error = simplify(factor * (tail - deep)).as_leading_term(eps)
    assert shallow_error.as_powers_dict()[eps] == target
    assert deep_error.as_powers_dict()[eps] >= target + 1
    assert simplify(factor * tail - exact) == 0


def test_tail_order_budget_adds_external_prefactor_poles():
    eps = symbols("eps")
    # A generic regular tail coefficient and a C_N with one pole, plus an
    # independent two-pole prefactor, require target+3 tail coefficients.
    factor = 1 / (eps**3 * (1 + eps))
    tail = (1 + 2 * eps) / (1 + 3 * eps)
    target = 2
    deep = tail.series(eps, 0, target + 3 + 1).removeO()
    shallow = tail.series(eps, 0, target + 3).removeO()
    assert simplify(factor * (tail - deep)).as_leading_term(eps).as_powers_dict()[eps] == 3
    assert simplify(factor * (tail - shallow)).as_leading_term(eps).as_powers_dict()[eps] == 2
