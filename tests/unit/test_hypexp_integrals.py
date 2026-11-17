import mpmath as mp
import pytest
import sympy as s

from hypra import HPL, UnsupportedExpressionError, beta_log_moment, hypergeometric_log_moment
from hypra.symbolic import reduce_known_hpl


@pytest.mark.parametrize(
    ("powers", "expected"),
    [
        ((3, 0, 2), -s.Rational(2, 27)),
        ((1, 1, 0), 2 - s.zeta(2)),
        ((0, 1, -1), -s.zeta(2)),
        ((0, 2, -2), 2 * s.zeta(2)),
    ],
)
def test_beta_log_moments(powers, expected):
    assert s.simplify(beta_log_moment(*powers) - expected) == 0


def test_moving_kernel_keeps_denominator_sign_and_symbolic_argument():
    z = s.Symbol("z")
    formulas = {1: s.log(1 - z) / z, 2: 1 / (1 - z), 3: (1 - (1 - z) ** -2) / (2 * z)}
    for q, expected in formulas.items():
        actual = reduce_known_hpl(hypergeometric_log_moment(0, 0, 0, 0, q, z))
        assert s.simplify(actual - expected) == 0


def test_moving_kernel_against_independent_quadrature():
    z = s.Rational(2, 5)
    actual = hypergeometric_log_moment(2, 1, 0, 2, 2, z)
    with mp.workdps(35):
        x = mp.mpf(2) / 5
        expected = mp.quad(
            lambda u: u**2 * mp.log(u) * mp.log1p(-x * u) ** 2 / (x * u - 1) ** 2, [0, 1]
        )
        assert actual.has(HPL)
        assert abs(mp.mpf(str(actual.evalf(25))) - expected) < mp.mpf("1e-24")


def test_endpoints_and_argument_independent_case():
    assert hypergeometric_log_moment(0, 1, 0, 0, 1, 0) == 1
    assert hypergeometric_log_moment(0, 1, 1, 1, 1, 0) == 0
    assert hypergeometric_log_moment(0, 1, 1, 0, 1, 1) == -s.zeta(3)
    assert hypergeometric_log_moment(2, 1, 0, 0, 0, 2) == -s.Rational(1, 9)


def test_unsupported_inputs():
    with pytest.raises(UnsupportedExpressionError, match="diverges at zero"):
        beta_log_moment(0, 1, -2)
    with pytest.raises(TypeError, match="integer"):
        beta_log_moment(True, 0, 0)
    with pytest.raises(UnsupportedExpressionError, match="diverges at one"):
        hypergeometric_log_moment(0, 0, 1, 1, 1, 1)
    with pytest.raises(UnsupportedExpressionError, match="0 < z < 1"):
        hypergeometric_log_moment(0, 0, 0, 1, 1, 2)


@pytest.mark.parametrize(
    "powers",
    [
        (0, 0, 1, 1, 0),
        (0, 0, 1, 0, 1),
        (2, 1, 1, 1, 2),
        (2, 1, 2, 1, 0),
        (3, 0, 1, 0, 4),
    ],
)
def test_mixed_moments_against_quadrature(powers):
    p, n, m, k, q = powers
    z = s.Symbol("z")
    exact = hypergeometric_log_moment(*powers, z)
    assert not exact.has(s.Integral, s.hyper, s.Float)
    with mp.workdps(40):
        x = mp.mpf(2) / 5
        expected = mp.quad(
            lambda u: (
                u**p * mp.log(u) ** n * mp.log1p(-u) ** m * mp.log1p(-x * u) ** k / (x * u - 1) ** q
            ),
            [0, 1],
        )
        actual = exact.subs(z, s.Rational(2, 5)).evalf(28)
        assert abs(mp.mpf(str(actual)) - expected) < mp.mpf("1e-26")


def test_mixed_moment_large_polynomial_degree():
    p = 36
    actual = beta_log_moment(0, 1, p)
    assert actual == -s.harmonic(p + 1) / (p + 1)


@pytest.mark.parametrize("n,p", [(8, 128), (12, 4096), (32, 0)])
def test_beta_origin_logs_do_not_require_hpl_endpoint_tables(n, p):
    # Direct differentiation of integral_0^1 u**p du = 1/(p+1).
    expected = s.S.NegativeOne**n * s.factorial(n) / s.Integer(p + 1) ** (n + 1)
    assert beta_log_moment(n, 0, p) == expected


@pytest.mark.parametrize("n", [7, 8, 11])
def test_beta_negative_power_high_weight_against_defining_log_series(n):
    # Expand log(1-u)=-sum_{k>=1} u**k/k, then integrate every log(u)**n term.
    expected = s.S.NegativeOne ** (n + 1) * s.factorial(n) * s.zeta(n + 2)
    assert s.simplify(beta_log_moment(n, 1, -1) - expected) == 0


def test_beta_large_degree_against_independent_gamma_derivatives():
    a, b = s.symbols("a b", positive=True)
    beta = s.gamma(a) * s.gamma(b) / s.gamma(a + b)
    expected = s.diff(beta, a, 3, b, 2).subs({a: 129, b: 1})
    assert s.simplify(beta_log_moment(3, 2, 128) - expected) == 0


@pytest.mark.parametrize("powers", [(3, 5, 0), (1, 8, -5), (0, 8, -8)])
def test_beta_high_weight_removable_poles_against_quadrature(powers):
    n, m, p = powers
    actual = beta_log_moment(*powers)
    assert not actual.has(HPL, s.Integral, s.Float, s.zoo, s.nan)
    with mp.workdps(60):
        expected = mp.quad(lambda u: u**p * mp.log(u) ** n * mp.log1p(-u) ** m, [0, 1])
        assert abs(mp.mpf(str(actual.evalf(45))) - expected) < mp.mpf("1e-40")


def test_beta_moment_budgets_and_convergence_remain_explicit():
    for powers in ((33, 0, 0), (0, 0, 4097)):
        with pytest.raises(UnsupportedExpressionError, match="budget"):
            beta_log_moment(*powers)
    with pytest.raises(UnsupportedExpressionError, match="diverges at zero"):
        beta_log_moment(0, 8, -9)


def test_negative_argument_reflection_preserves_all_log_factors():
    actual = hypergeometric_log_moment(1, 1, 1, 1, 1, -s.Rational(3, 2))
    with mp.workdps(40):
        x = -mp.mpf(3) / 2
        expected = mp.quad(
            lambda u: u * mp.log(u) * mp.log1p(-u) * mp.log1p(-x * u) / (x * u - 1),
            [0, 1],
        )
        assert abs(mp.mpf(str(actual.evalf(28))) - expected) < mp.mpf("1e-26")


@pytest.mark.parametrize("powers", [(0, 0, 1, 1, 0), (2, 1, 0, 1, 1)])
@pytest.mark.parametrize("point", [s.Rational(1, 2) + s.I / 3, 1 - s.I])
def test_complex_log_moments_follow_principal_radial_branch(powers, point):
    actual = hypergeometric_log_moment(*powers, point, argument_domain="principal-cut-plane")
    p, n, m, k, q = powers
    with mp.workdps(45):
        z = mp.mpc(str(s.re(point).evalf(45)), str(s.im(point).evalf(45)))
        expected = mp.quad(
            lambda u: (
                u**p
                * mp.log(u) ** n
                * mp.log1p(-u) ** m
                * mp.log(1 - z * u) ** k
                / (z * u - 1) ** q
            ),
            [0, 1],
        )
        value = actual.evalf(32)
        computed = mp.mpc(str(s.re(value)), str(s.im(value)))
        assert abs(computed - expected) < mp.mpf("1e-28")


def test_complex_log_moment_domain_is_explicit_and_keeps_endpoint_convergence():
    with pytest.raises(UnsupportedExpressionError, match="0 < z < 1"):
        hypergeometric_log_moment(0, 0, 0, 1, 0, s.I)
    with pytest.raises(UnsupportedExpressionError, match="exact finite numeric"):
        hypergeometric_log_moment(
            0, 0, 0, 1, 0, s.Symbol("z"), argument_domain="principal-cut-plane"
        )
    with pytest.raises(UnsupportedExpressionError, match="0 < z < 1"):
        hypergeometric_log_moment(0, 0, 0, 1, 0, 2, argument_domain="principal-cut-plane")
    with pytest.raises(UnsupportedExpressionError, match="diverges at one"):
        hypergeometric_log_moment(0, 0, 0, 0, 1, 1, argument_domain="principal-cut-plane")
    assert hypergeometric_log_moment(
        0, 1, 1, 0, 1, 1, argument_domain="principal-cut-plane"
    ) == -s.zeta(3)
    with pytest.raises(ValueError, match="argument_domain"):
        hypergeometric_log_moment(0, 0, 0, 1, 0, s.I, argument_domain="unknown")
