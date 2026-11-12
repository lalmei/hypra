"""Offline API contract tests. FakeBackend is not Wolfram validation."""

import pytest
import sympy as s

from hypra import (
    HPL,
    IncompleteExpansionError,
    UnsupportedExpressionError,
    encode_expr,
    expand_epsilon,
)
from hypra.planning import plan_expansion

E, W = s.symbols("eps w")
P = s.hyper([1, -E], [1 - 2 * E], W)


class FakeBackend:
    def __init__(self, polynomial, **overrides):
        self.polynomial, self.overrides, self.requests = polynomial, overrides, []

    def expand(self, request):
        self.requests.append(request)
        n = request["through"]
        return {
            "result": encode_expr(self.polynomial),
            "achieved_through": n,
            "remainder_order": n + 1,
            "complete": True,
            "backend": {"kind": "offline-test-double"},
            **self.overrides,
        }


def test_full_prefactor_requests_deeper_order():
    f = s.gamma(-E) ** 2 * P
    b = FakeBackend(E**-2 + 2 * E)
    result = expand_epsilon(f, E, through=1, backend=b)
    assert b.requests[0]["internal_order"] == 3
    assert result.coefficient(-2) == 1
    assert result.coefficient(1) == 2
    assert result.provenance["injected_backend"] is True
    assert result.provenance["plan"]["pole_bound"] == 2
    assert s.sympify(result) == E**-2 + 2 * E + s.Order(E**2)


def test_gamma_ratio_negative_powers_and_request_order():
    f = (1 - 2 * E) * s.gamma(1 - E) ** 2 / (E**2 * s.gamma(1 - 2 * E)) * P
    assert plan_expansion(f, E, 3).internal_order == 5


def test_lower_parameter_pole_budget():
    h = s.hyper([1, 1], [E], W, evaluate=False)
    p = plan_expansion(h * E**-2, E, 3)
    assert p.pole_bound == 3
    assert p.internal_order == 6
    assert p.valuation_lower_bound == -3


def test_product_of_hypergeometric_poles_is_conservative():
    h = s.hyper([1, 1], [E], W, evaluate=False)
    p = plan_expansion(h**2, E, 2)
    assert p.pole_bound == 2
    assert p.internal_order == 4


def test_negative_through():
    b = FakeBackend(E**-2)
    result = expand_epsilon(P / E**2, E, through=-1, backend=b)
    assert result.expression == E**-2 + s.Order(1, (E, 0))
    assert b.requests[0]["internal_order"] == 1
    with pytest.raises(ValueError):
        result.coefficient(0)


def test_known_zero_in_denominator_leading_coefficient_requires_assumption():
    f = P / (W + E)
    with pytest.raises(UnsupportedExpressionError, match="nonzero assumption"):
        plan_expansion(f, E, 3)
    assert plan_expansion(f, E, 3, s.Ne(W, 0)).internal_order == 3


def test_positive_domain_for_regulator_powers():
    f = W ** (-1 - 2 * E) * (1 - W) ** (-E) * P
    domain = s.And(W > 0, W < 1)
    assert plan_expansion(f, E, 3, domain).internal_order == 3
    with pytest.raises(UnsupportedExpressionError, match="positive"):
        plan_expansion(f, E, 3)


@pytest.mark.parametrize(
    "expr",
    [
        s.exp(1 / E),
        s.log(E),
        s.sqrt(E),
        1 / (1 + P),
        s.gamma(W + E),
        s.hyper([1, E], [1 + E], W + E),
    ],
)
def test_ambiguous_or_nonlaurent_cases_rejected(expr):
    with pytest.raises(UnsupportedExpressionError):
        expand_epsilon(expr, E, backend=FakeBackend(0))


def test_half_integer_input_is_admitted_not_proven_backend_supported():
    f = s.hyper([s.Rational(1, 2) + E, 1], [s.Rational(3, 2)], W)
    assert plan_expansion(f, E, 2).internal_order == 2


def test_incomplete_response_rejected():
    with pytest.raises(IncompleteExpansionError):
        expand_epsilon(P, E, backend=FakeBackend(1, complete=False))
    with pytest.raises(IncompleteExpansionError):
        expand_epsilon(P, E, backend=FakeBackend(1, achieved_through=2))


@pytest.mark.parametrize("out", [P, E**4, 1 / E, s.exp(E), s.log(E), s.Symbol("alien")])
def test_invalid_output_rejected(out):
    with pytest.raises(IncompleteExpansionError):
        expand_epsilon(P, E, backend=FakeBackend(out))


def test_exact_input_only():
    with pytest.raises(UnsupportedExpressionError):
        expand_epsilon(P * s.Float("0.1"), E, backend=FakeBackend(1))
    with pytest.raises(TypeError):
        expand_epsilon("gamma(eps)", E)
    with pytest.raises(TypeError):
        expand_epsilon(P, E, through=True, backend=FakeBackend(1))


def test_hpl_result_and_latex_roundtrip():
    b = FakeBackend(1 + E * HPL((0, 1), W))
    result = expand_epsilon(P, E, backend=b)
    assert result.coefficient(1) == HPL((0, 1), W)
    assert r"H_{0,1}" in result.latex()
    assert result.to_dict()["remainder_order"] == 4


def test_assumed_regulator_symbol_is_restored():
    eps = s.Symbol("eps", real=True)
    result = expand_epsilon(1 + eps, eps, backend=FakeBackend(1 + s.Symbol("eps")))
    assert result.coefficient(1) == 1
    assert result.polynomial.free_symbols == {eps}


def test_same_name_different_symbol_assumptions_rejected():
    with pytest.raises(UnsupportedExpressionError, match="ambiguous"):
        expand_epsilon(E + s.Symbol("eps", real=True), E, backend=FakeBackend(1))


@pytest.mark.parametrize("invalid", [{"head": "List", "args": []}, {"head": "True"}])
def test_nonscalar_result_is_explicit_incomplete_failure(invalid):
    with pytest.raises(IncompleteExpansionError, match="scalar"):
        expand_epsilon(P, E, backend=FakeBackend(1, result=invalid))
