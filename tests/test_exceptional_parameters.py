"""Safety regressions for exceptional parameter strata and limit order."""
import pytest
import sympy as s

from hypra import NativeConfig, expand_epsilon
from hypra.endpoints import rewrite_endpoints
from hypra.errors import UnsupportedExpressionError
from hypra.planning import plan_expansion

e,z=s.symbols('e z', real=True)
domain=s.And(z>0,z<1)
config=NativeConfig(use_cache=False)


@pytest.mark.parametrize('argument',[0,z,1])
def test_fixed_singular_equal_pair_never_canceled(argument):
    h=s.hyper([-1,2+e],[-1],argument,evaluate=False)
    with pytest.raises(UnsupportedExpressionError,match='lower'):
        rewrite_endpoints(h,e)
    with pytest.raises(UnsupportedExpressionError,match='lower'):
        expand_epsilon(h,e,through=1,assumptions=domain,config=config)


def test_regulator_independent_singularity_rejected_by_planner():
    h=s.hyper([2,3],[-1],z,evaluate=False)
    with pytest.raises(UnsupportedExpressionError,match='lower'):
        plan_expansion(1+h,e,1,domain)


def test_equal_moving_singular_pairs_cancel_on_punctured_neighborhood():
    h=s.hyper([-1+e,2+e],[-1+e],z,evaluate=False)
    reduced=rewrite_endpoints(h,e)
    assert reduced.ap==(2+e,) and reduced.bq==()
    assert plan_expansion(reduced,e,2,domain).pole_bound==0
    actual=expand_epsilon(h,e,through=2,assumptions=domain,config=config).polynomial
    expected=(1-e*s.log(1-z)+e**2*s.log(1-z)**2/2)/(1-z)**2
    assert s.simplify(actual-expected)==0


def test_termination_with_moving_singular_lower_has_real_pole():
    h=s.hyper([-2,1+e],[-e],z,evaluate=False)
    finite=1+2*(1+e)*z/e-(1+e)*(2+e)*z*z/(e*(1-e))
    rewritten=rewrite_endpoints(h,e)
    assert s.cancel(rewritten-finite)==0
    result=expand_epsilon(h,e,through=1,assumptions=domain,config=config).polynomial
    expected=s.series(finite,e,0,2).removeO()
    assert s.simplify(result-expected)==0


def test_zero_upper_still_requires_valid_lower_but_moving_lower_is_valid():
    good=s.hyper([0,2],[-e],z,evaluate=False)
    assert rewrite_endpoints(good,e)==1
    bad=s.hyper([0,2],[-1],z,evaluate=False)
    with pytest.raises(UnsupportedExpressionError):
        rewrite_endpoints(bad,e)


def test_equal_pairs_preserve_convergent_gauss_zero():
    h=s.hyper([2+e,-s.Rational(1,2)],[2+e],1,evaluate=False)
    assert rewrite_endpoints(h,e)==0


def test_cancellation_reduces_rank_at_unit_endpoint():
    h=s.hyper([2+e,1,1],[2+e,3],1,evaluate=False)
    assert rewrite_endpoints(h,e)==2


def test_degenerate_contiguous_pair_limit_recovers_exact_binomial():
    from hypra.canonical import build_canonical
    from hypra.contiguous import expand_contiguous, plan_contiguous
    from hypra.symbolic import reduce_known_hpl
    upper,lower=[(s.Integer(2),s.S.Zero),(s.S.One,s.S.Zero)],[(s.Integer(2),s.S.Zero)]
    plan=plan_contiguous(upper,lower,[s.S.One,s.S.One],[s.S.One])
    assert plan.perturbation_used
    result=expand_contiguous(upper,lower,[s.S.One,s.S.One],[s.S.One],z,e,0,build_canonical)
    assert s.simplify(reduce_known_hpl(result)-1/(1-z))==0
    assert result.free_symbols <= {z}
