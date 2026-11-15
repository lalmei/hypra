import pytest
import sympy as s

from hypra import NativeConfig, expand_epsilon
from hypra.errors import UnsupportedExpressionError
from hypra.native.endpoints import rewrite_endpoints

e = s.Symbol("epsilon", real=True)


def test_general_integer_endpoint_constant_and_regulator_coefficient():
    # C_n=(1+e)/[(n+1)(n+1+e)], hence sum gives
    # zeta(2)+e*(zeta(2)-zeta(3))+e²*(zeta(4)-zeta(3)).
    expr = s.hyper([1, 1, 1 + e], [2, 2 + e], 1)
    result = expand_epsilon(expr, e, through=2, config=NativeConfig(use_cache=False)).polynomial
    expected = s.zeta(2) + e * (s.zeta(2) - s.zeta(3)) + e**2 * (s.zeta(4) - s.zeta(3))
    assert s.simplify(result - expected) == 0
    constant = expand_epsilon(
        s.hyper([1, 1, 1], [2, 2], 1), e, through=0, config=NativeConfig(use_cache=False)
    ).polynomial
    assert constant == s.zeta(2)


def test_endpoint_prefactor_budget_and_nonpositive_excess_rejection():
    expr = s.hyper([1, 1, 1 + e], [2, 2 + e], 1) / e
    result = expand_epsilon(expr, e, through=0, config=NativeConfig(use_cache=False)).polynomial
    assert s.simplify(result - (s.zeta(2) / e + s.zeta(2) - s.zeta(3))) == 0
    for h in (
        s.hyper([1, 1, 2], [2, 2], 1, evaluate=False),
        s.hyper([1, 1, 3 + e], [2, 2], 1, evaluate=False),
    ):
        with pytest.raises(UnsupportedExpressionError):
            rewrite_endpoints(h, e)


def test_moving_lower_pole_survives_uniform_endpoint_limit():
    # Isolate n=1, the only coefficient with an epsilon pole. The n>=2
    # constant term is a convergent rational telescoping sum.
    n = s.Symbol("n", integer=True, positive=True)
    tail = s.summation(6 / ((n - 1) * (n + 1) * (n + 2) * (n + 3)), (n, 2, s.oo))
    expected = 1 + s.Rational(1, 4) - tail - 1 / (4 * e)
    result = expand_epsilon(
        s.hyper([-1 + e, 1, 1], [e, 4], 1), e, through=0, config=NativeConfig(use_cache=False)
    ).polynomial
    assert s.simplify(result - expected) == 0


def test_public_endpoint_rewrites_before_pole_planning():
    eps = s.Symbol("eps", real=True)
    result = expand_epsilon(s.hyper([1, 1], [2 - eps], 1), eps, through=3, assumptions=eps < 0)
    assert s.expand(result.polynomial - (1 - 1 / eps)) == 0
    assert result.provenance["plan"]["pole_bound"] == 1
    assert result.provenance["input"] != result.provenance["planned_input"]
