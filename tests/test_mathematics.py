"""Offline mathematical reference tests, not real-Wolfram backend tests."""

import importlib.util
from pathlib import Path

import sympy as s

MODULE = Path(__file__).resolve().parents[1] / "scripts" / "validate_integral.py"
SPEC = importlib.util.spec_from_file_location("validate_integral", MODULE)
validation = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validation)


def test_p_coefficient_derivation_against_finite_pochhammer_ratios():
    eps = s.Symbol("epsilon")
    for n in range(1, 10):
        ratio = s.rf(-eps, n) / s.rf(1 - 2 * eps, n)
        regular = s.cancel(ratio / eps)
        assert regular.subs(eps, 0) == -s.Rational(1, n)
        expected_second = -(s.harmonic(n - 1) + s.Rational(2, n)) / n
        assert s.simplify(s.diff(regular, eps).subs(eps, 0) - expected_second) == 0


def test_q_ratio_telescopes_exactly():
    eps = s.Symbol("epsilon")
    for n in range(1, 10):
        ratio = s.rf(-2 * eps, n) / s.rf(1 - 2 * eps, n)
        assert s.cancel(ratio + 2 * eps / (n - 2 * eps)) == 0


def test_beta_derivative_and_reflection_moments():
    moments = validation.derive_moments()
    assert s.simplify(moments["log(w)*log(1-w)^2"] + s.zeta(4) / 2) == 0
    assert s.simplify(moments["log(1-w)*Li2(1-w)"] + 3 * s.zeta(4) / 4) == 0


def test_j4_derived_not_imported_from_nnlo_table():
    derived = validation.derive_j4()
    assert derived["coefficient_2"] == 2 * s.zeta(3)
    assert s.simplify(derived["coefficient_3"] - 11 * s.pi**4 / 60) == 0
    assert derived["J4"].getO() == s.Order(validation.EPS**4, validation.EPS)


def _rational_taylor_coefficients(expression, variable, order):
    """Exact polynomial division, independent of the HPL differential equation."""
    numerator, denominator = s.fraction(s.cancel(expression))
    numerator, denominator = s.Poly(numerator, variable), s.Poly(denominator, variable)
    d0 = denominator.nth(0)
    assert d0 != 0
    result = []
    for k in range(order + 1):
        previous = sum(denominator.nth(j) * result[k - j] for j in range(1, k + 1))
        result.append(s.cancel((numerator.nth(k) - previous) / d0))
    return result


def test_generic_first_order_ode_recursion_matches_exact_pochhammer_series():
    """Test regulator orders 1..5 and argument powers 1..6 for four families."""
    eps = s.Symbol("epsilon")
    max_order, max_power = 5, 6
    parameters = [(-1, -2), (-2, -2), (s.Rational(2, 3), -s.Rational(1, 4)), (0, 3)]
    for b, c in parameters:
        # Integrating -c*f/z+(b-c)*f/(1-z) coefficientwise is the HPL
        # prepend-0/prepend-1 recursion, with the regular zero boundary value.
        coefficients = {1: [s.S.Zero] + [s.Rational(b, n) for n in range(1, max_power + 1)]}
        for k in range(2, max_order + 1):
            previous = coefficients[k - 1]
            coefficients[k] = [s.S.Zero] + [
                s.cancel((-c * previous[n] + (b - c) * sum(previous[1:n])) / n)
                for n in range(1, max_power + 1)
            ]
        for n in range(1, max_power + 1):
            exact = _rational_taylor_coefficients(
                s.rf(b * eps, n) / s.rf(1 + c * eps, n), eps, max_order
            )
            assert exact[0] == 0
            for k in range(1, max_order + 1):
                assert s.simplify(coefficients[k][n] - exact[k]) == 0


def test_half_integer_even_binomial_identity():
    a = s.Symbol("a")
    for n in range(8):
        even_binomial = s.rf(2 * a, 2 * n) / s.factorial(2 * n)
        hyper_coefficient = (
            s.rf(a, n)
            * s.rf(a + s.Rational(1, 2), n)
            / (s.rf(s.Rational(1, 2), n) * s.factorial(n))
        )
        assert s.cancel(even_binomial - hyper_coefficient) == 0


def test_half_integer_logarithmic_coefficients_match_hypergeometric_series():
    eps, t = s.symbols("epsilon t")
    a = s.Rational(2, 3)
    # Use t=sqrt(z). Odd terms cancel; coefficients through t^8 suffice.
    minus_log = -sum(t**j / j for j in range(1, 9))
    plus_log = sum((-1) ** (j + 1) * t**j / j for j in range(1, 9))
    for k in range(1, 4):
        f_k = s.expand((-2 * a) ** k * (minus_log**k + plus_log**k) / (2 * s.factorial(k)))
        for n in range(1, 5):
            exact = _rational_taylor_coefficients(
                s.rf(a * eps, n)
                * s.rf(s.Rational(1, 2) + a * eps, n)
                / (s.rf(s.Rational(1, 2), n) * s.factorial(n)),
                eps,
                k,
            )
            assert s.simplify(f_k.coeff(t, 2 * n) - exact[k]) == 0
        assert all(f_k.coeff(t, j) == 0 for j in range(1, 9, 2))


def test_half_integer_telescoping_coefficients():
    eps = s.Symbol("epsilon")
    b = s.Rational(3, 4)
    a0 = s.Rational(1, 2)
    for n in range(7):
        exact_ratio = s.rf(a0 + b * eps, n) / s.rf(1 + a0 + b * eps, n)
        assert s.cancel(exact_ratio - (a0 + b * eps) / (n + a0 + b * eps)) == 0
        exact = _rational_taylor_coefficients(exact_ratio, eps, 5)
        assert exact[0] == a0 / (n + a0)
        for k in range(1, 6):
            odd_polylog_coefficient = (
                (-1) ** (k - 1) * b**k * ((n + a0) ** (-k) - a0 * (n + a0) ** (-k - 1))
            )
            assert s.cancel(exact[k] - odd_polylog_coefficient) == 0
