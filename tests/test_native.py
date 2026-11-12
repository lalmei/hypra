"""Real native-Python calculations, with independent coefficient/numeric checks."""

import json

import mpmath as mp
import pytest
import sympy as s

from hypra import NativeConfig, UnsupportedExpressionError, backend_info, expand_epsilon
from hypra.native import expand_atom
from hypra.serialization import encode_expr
from hypra.symbolic import reduce_known_hpl

E, W = s.symbols("eps w")
DOMAIN = s.And(W > 0, W < 1)
NO_CACHE = NativeConfig(use_cache=False)


def test_native_default_needs_no_external_executable():
    assert backend_info()["available"] is True
    assert backend_info()["backend"] == "native-sympy"


def test_p_q_acceptance_computed_natively():
    p = expand_epsilon(s.hyper([1, -E], [1 - 2 * E], W), E, 2, DOMAIN, config=NO_CACHE)
    q = expand_epsilon(s.hyper([1, -2 * E], [1 - 2 * E], 1 - W), E, 3, DOMAIN, config=NO_CACHE)
    assert s.simplify(p.coefficient(1) - s.log(1 - W)) == 0
    assert s.simplify(p.coefficient(2) + 2 * s.polylog(2, W) + s.log(1 - W) ** 2 / 2) == 0
    assert s.simplify(q.coefficient(1) - 2 * s.log(W)) == 0
    assert s.simplify(q.coefficient(2) + 4 * s.polylog(2, 1 - W)) == 0
    assert s.simplify(q.coefficient(3) + 8 * s.polylog(3, 1 - W)) == 0
    assert p.provenance["backend"]["backend"] == "native-sympy"
    assert p.provenance["injected_backend"] is False


def test_generic_recurrence_satisfies_differential_equation():
    b, c = s.Rational(3, 2), s.Rational(-1, 2)
    h = s.hyper([1, b * E], [1 + c * E], W)
    polynomial = expand_atom(h, E, 4, DOMAIN)
    residual = s.expand(W * (1 - W) * s.diff(polynomial, W) + E * (c - b * W) * polynomial - c * E)
    for k in range(5):
        assert s.simplify(s.expand_func(reduce_known_hpl(residual.coeff(E, k)))) == 0


def test_q_all_orders_not_example_lookup():
    b = s.Rational(3, 2)
    h = s.hyper([1, b * E], [1 + b * E], W)
    result = expand_epsilon(h, E, 6, DOMAIN, config=NO_CACHE)
    for k in range(1, 7):
        expected = (-1) ** (k - 1) * b**k * s.polylog(k, W)
        assert s.simplify(s.expand_func(result.coefficient(k) - expected)) == 0


def test_pole_prefactor_uses_q_cubic_term():
    q = s.hyper([1, -2 * E], [1 - 2 * E], 1 - W)
    result = expand_epsilon(s.gamma(E) * q, E, 2, DOMAIN, config=NO_CACHE)
    expected_q = (
        1 + 2 * E * s.log(W) - 4 * E**2 * s.polylog(2, 1 - W) - 8 * E**3 * s.polylog(3, 1 - W)
    )
    expected = s.series(s.gamma(E) * expected_q, E, 0, 3).removeO().expand()
    assert s.simplify(result.polynomial - expected) == 0
    assert result.provenance["plan"]["internal_order"] == 3


def test_half_integer_telescoping_family():
    h = s.hyper([1, s.Rational(1, 2) + E], [s.Rational(3, 2) + E], W)
    result = expand_epsilon(h, E, 2, DOMAIN, config=NO_CACHE)
    assert s.simplify(result.coefficient(0).subs(W, s.Rational(1, 4)) - s.log(3)) == 0
    # Direct rational coefficients of z^n verify epsilon coefficients independently.
    for n in range(1, 5):
        a = s.Rational(1, 2) + E
        expected = s.series(a / (a + n), E, 0, 3).removeO()
        actual = s.series(result.polynomial, W, 0, n + 1).removeO().expand().coeff(W, n)
        assert s.simplify(actual - expected) == 0


def test_half_integer_even_binomial_family():
    h = s.hyper([E, s.Rational(1, 2) + E], [s.Rational(1, 2)], W)
    result = expand_epsilon(h, E, 3, DOMAIN, config=NO_CACHE)
    t = s.sqrt(W)
    exact = ((1 - t) ** (-2 * E) + (1 + t) ** (-2 * E)) / 2
    expected = s.series(exact, E, 0, 4).removeO()
    assert s.simplify(result.polynomial - expected) == 0


def test_native_domain_and_unsupported_families_fail():
    with pytest.raises(UnsupportedExpressionError, match="0 <"):
        expand_epsilon(s.hyper([1, -E], [1 - 2 * E], W), E, config=NO_CACHE)
    with pytest.raises(UnsupportedExpressionError, match="family"):
        expand_epsilon(s.hyper([2, E], [1 - E], W), E, assumptions=DOMAIN, config=NO_CACHE)
    with pytest.raises(UnsupportedExpressionError, match="limit"):
        expand_epsilon(s.hyper([1, E], [1 - E], W), E, through=13, assumptions=DOMAIN)


def test_cache_uses_exact_inputs_and_orders(tmp_path):
    config = NativeConfig(cache_dir=str(tmp_path))
    f = s.hyper([1, -E], [1 - 2 * E], W)
    first = expand_epsilon(f, E, 2, DOMAIN, config=config)
    second = expand_epsilon(f, E, 2, DOMAIN, config=config)
    assert first.polynomial == second.polynomial
    assert first.provenance["cache_hit"] is False
    assert second.provenance["cache_hit"] is True
    assert len(list(tmp_path.glob("*.json"))) == 1
    expand_epsilon(f, E, 3, DOMAIN, config=config)
    assert len(list(tmp_path.glob("*.json"))) == 2
    path = next(tmp_path.glob("*.json"))
    data = json.loads(path.read_text())
    data["response"]["result"] = encode_expr(s.Integer(123))
    path.write_text(json.dumps(data))
    # Mismatched checksum is ignored, never trusted as a coefficient result.
    fresh = expand_epsilon(f, E, 2, DOMAIN, config=config)
    assert fresh.polynomial != 123


def test_secondary_numerical_q_remainder():
    q = s.hyper([1, -2 * E], [1 - 2 * E], W)
    result = expand_epsilon(q, E, 3, DOMAIN, config=NO_CACHE)
    fn = s.lambdify((E, W), result.polynomial, "mpmath")
    with mp.workdps(40):
        errors = []
        for eps in (mp.mpf("-.01"), mp.mpf("-.005")):
            exact = mp.hyp2f1(1, -2 * eps, 1 - 2 * eps, mp.mpf(".3"))
            errors.append(abs(exact - fn(eps, mp.mpf(".3"))))
        assert 15 < errors[0] / errors[1] < 17


@pytest.mark.parametrize(
    "expr,through,expected",
    [
        (E, -1, s.S.Zero),
        (E, -3, s.S.Zero),
        (E**-2, -2, E**-2),
        (E**-2 + E, -1, E**-2),
        (s.gamma(E), -1, E**-1),
    ],
)
def test_negative_requested_orders(expr, through, expected):
    assert expand_epsilon(expr, E, through, config=NO_CACHE).polynomial == expected


def test_contradictory_domain_is_rejected_even_without_hyper():
    with pytest.raises(UnsupportedExpressionError, match="Contradictory"):
        expand_epsilon(1 + E, E, assumptions=s.And(W > 0, W < 0, evaluate=False), config=NO_CACHE)


def test_actual_native_j4_integrand_before_proved_integration():
    p = s.hyper([1, -E], [1 - 2 * E], W)
    q = s.hyper([1, -2 * E], [1 - 2 * E], 1 - W)
    integrand = W ** (-1 - 2 * E) * (1 - W) ** (-E) * (p - 1) * (q - 1)
    result = expand_epsilon(integrand, E, 3, DOMAIN, config=NO_CACHE)
    l, m = s.log(W), s.log(1 - W)
    assert s.simplify(result.coefficient(2) - 2 * l * m / W) == 0
    expected = (
        -4 * l * s.polylog(2, W) - 3 * l * m**2 - 4 * m * s.polylog(2, 1 - W) - 4 * l**2 * m
    ) / W
    assert s.simplify(result.coefficient(3) - expected) == 0
