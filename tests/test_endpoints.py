import pytest
import sympy as s

from hypra.endpoints import rewrite_endpoints
from hypra.errors import UnsupportedExpressionError

e, z = s.symbols('e z', real=True)


def test_terminating_exact_finite_sum_any_rank():
    expr = s.hyper([-3, 1+e, 2-e], [2+e], z)
    expected = sum(z**n*s.rf(-3,n)*s.rf(1+e,n)*s.rf(2-e,n)/(s.rf(2+e,n)*s.factorial(n)) for n in range(4))
    assert s.expand(rewrite_endpoints(expr,e)-expected) == 0
    assert s.simplify(rewrite_endpoints(expr.subs(z,1),e)-expected.subs(z,1)) == 0


def test_termination_keeps_moving_lower_poles():
    result = rewrite_endpoints(s.hyper([-1,2+e],[-e],z),e)
    assert s.simplify(result-(1+z*(2+e)/e)) == 0


def test_singular_fixed_lower_rejected_even_if_termination_precedes_it():
    with pytest.raises(UnsupportedExpressionError, match='lower'):
        rewrite_endpoints(s.hyper([-1],[-3],z),e)


def test_gauss_positive_excess_neighborhood():
    a,b,c = 1+e,2-e,4+3*e
    actual = rewrite_endpoints(s.hyper([a,b],[c],1),e)
    expected = s.gamma(c)*s.gamma(c-a-b)/(s.gamma(c-a)*s.gamma(c-b))
    assert s.simplify(actual-expected) == 0


def test_gauss_boundary_requires_proven_side_and_reveals_pole():
    expr = s.hyper([1,1],[2-e],1)
    with pytest.raises(UnsupportedExpressionError, match='convergence'):
        rewrite_endpoints(expr,e)
    actual = rewrite_endpoints(expr,e,e<0)
    assert s.simplify(actual-(e-1)/e) == 0
    with pytest.raises(UnsupportedExpressionError, match='convergence'):
        rewrite_endpoints(expr,e,e>0)


@pytest.mark.parametrize('expr', [s.hyper([1,1],[2],1), s.hyper([2,2],[3+e],1),
                                      s.hyper([2,2,2],[3,3],1)])
def test_no_analytic_continuation_of_divergent_or_unsupported_endpoint(expr):
    with pytest.raises(UnsupportedExpressionError):
        rewrite_endpoints(expr,e)


def test_regular_argument_unchanged_and_only_fixed_upper_terminates():
    expr = s.hyper([-1+e,2],[3],z)
    assert rewrite_endpoints(expr,e) == expr


def test_fixed_other_parameter_is_preserved():
    a=s.Symbol('a')
    expr=s.hyper([-2,a],[3],z)
    assert s.expand(rewrite_endpoints(expr,e)-(1-2*a*z/3+a*(a+1)*z*z/12)) == 0


def test_invalid_and_nonlinear_endpoint():
    with pytest.raises(UnsupportedExpressionError):
        rewrite_endpoints(s.hyper([1+e**2,1],[3],1),e)
    with pytest.raises(UnsupportedExpressionError):
        rewrite_endpoints(s.hyper([1,1],[3],1),e,s.false)


def test_public_endpoint_rewrites_before_pole_planning():
    from hypra import expand_epsilon
    eps = s.Symbol("eps", real=True)
    result = expand_epsilon(s.hyper([1, 1], [2-eps], 1), eps, through=3, assumptions=eps < 0)
    assert s.expand(result.polynomial - (1 - 1/eps)) == 0
    assert result.provenance["plan"]["pole_bound"] == 1
    assert result.provenance["input"] != result.provenance["planned_input"]
