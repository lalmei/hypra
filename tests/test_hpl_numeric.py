"""Numerical checks are independent approximations, not symbolic proof fixtures."""

import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import shuffle_words
from hypra.hpl_numeric import HPLNumericalConvergenceError, evaluate_hpl, hpl_evalf


def close_to(actual, expected, digits=45):
    with mp.workdps(digits + 20):
        assert abs(actual - expected) <= mp.power(10, -digits) * abs(expected)


@pytest.mark.parametrize("argument", [s.Rational(1, 5), "0.73", s.Rational(99, 100)])
def test_elementary_and_classical_words(argument):
    with mp.workdps(75):
        z = mp.mpf(str(s.N(argument, 75))) if isinstance(argument, s.Expr) else mp.mpf(argument)
        for word, expected in [
            ((), mp.mpf(1)),
            ((0,), mp.log(z)),
            ((0, 0, 0), mp.log(z) ** 3 / 6),
            ((1, 1), mp.log1p(-z) ** 2 / 2),
            ((-1, -1, -1), mp.log1p(z) ** 3 / 6),
            ((0, 1), mp.polylog(2, z)),
            ((0, 0, -1), -mp.polylog(3, -z)),
        ]:
            result = evaluate_hpl(word, argument, dps=55)
            close_to(result.value, expected, 50)
            assert result.tail_bound == 0
            assert result.terms == 0


@pytest.mark.parametrize("argument", [s.Rational(1, 10), s.Rational(3, 5)])
def test_logarithmic_trailing_zero_words(argument):
    with mp.workdps(80):
        z = mp.mpf(int(argument.p)) / int(argument.q)
        expected = {
            (1, 0): -mp.log(z) * mp.log1p(-z) - mp.polylog(2, z),
            (0, 1, 0): mp.log(z) * mp.polylog(2, z) - 2 * mp.polylog(3, z),
            (-1, 0): mp.log(z) * mp.log1p(z) + mp.polylog(2, -z),
            (1, 0, 0): -(mp.log(z) ** 2) * mp.log1p(-z) / 2
            - mp.log(z) * mp.polylog(2, z)
            + mp.polylog(3, z),
        }
        for word, value in expected.items():
            result = evaluate_hpl(word, argument, dps=50)
            close_to(result.value, value, 45)
            assert result.method == "logarithmic-power-series"
            assert result.terms > 0
            assert result.tail_bound <= mp.mpf("1e-58") * abs(result.value)
            # The independent identity's error is smaller than the truncation
            # majorant plus a separate (noncertified) roundoff allowance.
            assert abs(result.value - value) <= result.tail_bound + mp.mpf("1e-70")


def test_mixed_words_against_direct_single_integrals():
    with mp.workdps(65):
        z = mp.mpf("0.4")
        a = mp.quad(lambda t: mp.log1p(t) / (1 - t), [0, z])
        b = mp.quad(lambda t: -mp.log1p(-t) / (1 + t), [0, z])
        close_to(hpl_evalf((1, -1), "0.4", dps=50), a)
        close_to(hpl_evalf((-1, 1), "0.4", dps=50), b)


def test_arbitrary_word_shuffle_consistency():
    with mp.workdps(65):
        argument = s.Rational(2, 5)
        left, right = (0, -1), (1, 0)
        product = hpl_evalf(left, argument, dps=50) * hpl_evalf(right, argument, dps=50)
        shuffled = sum(
            count * hpl_evalf(word, argument, dps=50)
            for word, count in shuffle_words(left, right).items()
        )
        close_to(shuffled, product)


def test_high_precision_and_global_context_unchanged():
    old_precision = mp.mp.dps
    value = hpl_evalf((1, 0), s.Rational(1, 3), dps=85)
    assert mp.mp.dps == old_precision
    with mp.workdps(110):
        z = mp.mpf(1) / 3
        expected = -mp.log(z) * mp.log1p(-z) - mp.polylog(2, z)
        close_to(value, expected, 80)


def test_term_budget_failure_returns_no_partial_result():
    with pytest.raises(HPLNumericalConvergenceError, match="after 10 terms"):
        hpl_evalf((1, 0), "0.99", dps=50, max_terms=10)


@pytest.mark.parametrize("argument", [0, 1, -1, 2, float("nan"), float("inf"), s.I, 0.5 + 0j])
def test_endpoints_nonreal_and_outside_domain_rejected(argument):
    with pytest.raises((UnsupportedExpressionError, TypeError)):
        hpl_evalf((0, 1), argument)


@pytest.mark.parametrize("argument", [s.Symbol("x"), "__import__('os')", "1/3", object()])
def test_ambiguous_or_executable_input_rejected(argument):
    with pytest.raises((UnsupportedExpressionError, TypeError, ValueError)):
        hpl_evalf((0, 1), argument)


@pytest.mark.parametrize(
    "options", [{"dps": True}, {"dps": 1}, {"max_terms": 0}, {"max_terms": 1.5}]
)
def test_bad_config_rejected(options):
    with pytest.raises(ValueError):
        hpl_evalf((1, 0), "0.5", **options)


def test_invalid_word_rejected():
    with pytest.raises(UnsupportedExpressionError):
        hpl_evalf((0, 2), "0.5")
    with pytest.raises(UnsupportedExpressionError):
        hpl_evalf((1.0,), "0.5")


def test_public_numeric_exports():
    from hypra import HPLNumericResult
    from hypra import evaluate_hpl as public_evaluate

    result = public_evaluate((0,), "0.5", dps=20)
    assert isinstance(result, HPLNumericResult)
    assert result.method == "logarithm"


@pytest.mark.parametrize("argument", [s.I / 3, -s.I / 3, s.Rational(-1, 4) + s.I / 5])
def test_explicit_complex_disk_against_classical_identities(argument):
    with mp.workdps(80):
        a, b = argument.as_real_imag()
        z = mp.mpc(str(s.N(a, 80)), str(s.N(b, 80)))
        expected = {
            (0,): mp.log(z),
            (1, 1): mp.log1p(-z) ** 2 / 2,
            (-1, -1): mp.log1p(z) ** 2 / 2,
            (0, 0, 1): mp.polylog(3, z),
            (0, -1): -mp.polylog(2, -z),
            (1, 0): -mp.log(z) * mp.log1p(-z) - mp.polylog(2, z),
            (0, 1, 0): mp.log(z) * mp.polylog(2, z) - 2 * mp.polylog(3, z),
            (-1, 0): mp.log(z) * mp.log1p(z) + mp.polylog(2, -z),
        }
        for word, value in expected.items():
            result = evaluate_hpl(word, argument, domain="principal-disk", dps=55)
            close_to(result.value, value, 50)
            assert result.tail_bound >= 0
            assert "principal slit disk" in result.domain
            assert abs(result.value - value) <= result.tail_bound + mp.mpf("1e-75")


def test_complex_mixed_integral_and_shuffle():
    with mp.workdps(70):
        argument = s.Rational(1, 4) + s.I / 5
        z = mp.mpc(mp.mpf(1) / 4, mp.mpf(1) / 5)
        # Integrate along the straight radial path; no branch cut is crossed.
        expected = mp.quad(lambda t: z * mp.log1p(z * t) / (1 - z * t), [0, 1])
        close_to(hpl_evalf((1, -1), argument, domain="principal-disk", dps=55), expected, 50)
        left, right = (0, -1), (1, 0)
        def value(word):
            return hpl_evalf(word, argument, domain="principal-disk", dps=55)
        close_to(sum(count * value(word) for word, count in shuffle_words(left, right).items()),
                 value(left) * value(right), 50)


@pytest.mark.parametrize("argument", [0, 1, s.I, -s.Rational(1, 2), -0.5 + 0j,
                                       1 + 1j, complex(float("nan"), 0)])
def test_complex_disk_endpoints_boundary_cut_and_nonfinite_rejected(argument):
    with pytest.raises((UnsupportedExpressionError, ValueError)):
        hpl_evalf((0, 1), argument, domain="principal-disk")


def test_complex_stored_inputs_and_explicit_domain():
    with mp.workdps(70):
        for point in [0.25 + 0.125j, mp.mpc("0.25", "0.125")]:
            close_to(hpl_evalf((0,), point, domain="principal-disk", dps=50), mp.log(point))
        with pytest.raises(HPLNumericalConvergenceError):
            hpl_evalf((1, 0), s.I * s.Rational(99, 100), domain="principal-disk", max_terms=10)
    with pytest.raises(ValueError, match="domain"):
        hpl_evalf((0,), "0.5", domain="automatic")
