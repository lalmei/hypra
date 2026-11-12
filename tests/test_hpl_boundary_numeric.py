import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_boundary_numeric import evaluate_hpl_boundary


@pytest.mark.parametrize('point,region', [(s.Rational(-1,2),'-1<x<0'),(2,'x>1'),(-2,'x<-1')])
@pytest.mark.parametrize('side,sign', [('upper',1),('lower',-1)])
def test_logarithm_boundary_signs_exactly(point,region,side,sign):
    with mp.workdps(65):
        x = mp.mpf(int(point.p))/int(point.q) if isinstance(point,s.Rational) else mp.mpf(point)
        expected = {
            (0,): mp.log(abs(x)) + (sign*mp.pi*1j if x<0 else 0),
            (1,): -mp.log(abs(1-x)) + (sign*mp.pi*1j if x>1 else 0),
            (-1,): mp.log(abs(1+x)) + (sign*mp.pi*1j if x<-1 else 0),
        }
        for word,value in expected.items():
            result = evaluate_hpl_boundary(word,point,region=region,boundary=side,dps=45)
            assert abs(result.value-value) < mp.mpf('1e-45')
            assert result.boundary == side


def test_dilogarithm_positive_cut_and_conjugacy_without_offset():
    with mp.workdps(65):
        # Li2(2±i0)=pi²/4 ± i*pi*log(2), derived by real inversion.
        expected = mp.pi**2/4 + 1j*mp.pi*mp.log(2)
        upper = evaluate_hpl_boundary((0,1),2,region='x>1',boundary='upper',dps=45)
        lower = evaluate_hpl_boundary((0,1),2,region='x>1',boundary='lower',dps=45)
        assert abs(upper.value-expected) < mp.mpf('1e-45')
        assert abs(lower.value-mp.conj(upper.value)) < mp.mpf('1e-45')


def test_mixed_word_negative_cut_shuffle_and_component_diagnostics():
    with mp.workdps(60):
        def value(word):
            return evaluate_hpl_boundary(word,-2,region='x<-1',boundary='upper',dps=40)
        left,right = value((0,-1)),value((-1,0))
        product = value((0,)).value*value((-1,)).value
        assert abs(left.value+right.value-product) < mp.mpf('1e-40')
        assert all(tail>=0 for tail in left.component_tail_bounds)


@pytest.mark.parametrize('point,region', [(1,'x>1'),(-1,'x<-1'),(0,'-1<x<0'),
                                         (2,'x<-1'),(0.5,'0<x<1'),(s.I,'x>1')])
def test_invalid_exact_domain_rejected(point,region):
    with pytest.raises(UnsupportedExpressionError):
        evaluate_hpl_boundary((0,),point,region=region,boundary='upper')


def test_explicit_side_and_weight_limits():
    with pytest.raises(UnsupportedExpressionError):
        evaluate_hpl_boundary((0,),2,region='x>1',boundary='automatic')
    with pytest.raises(UnsupportedExpressionError):
        evaluate_hpl_boundary((0,)*7,2,region='x>1',boundary='upper')


def test_general_mixed_weight_three_conjugate_boundaries():
    with mp.workdps(55):
        upper = evaluate_hpl_boundary((-1,0,1),2,region='x>1',boundary='upper',dps=35)
        lower = evaluate_hpl_boundary((-1,0,1),2,region='x>1',boundary='lower',dps=35)
        assert abs(upper.value-mp.conj(lower.value)) < mp.mpf('1e-35')


def test_numerical_failure_and_global_context():
    from hypra.hpl_numeric import HPLNumericalConvergenceError
    original = mp.mp.dps
    with pytest.raises(HPLNumericalConvergenceError):
        evaluate_hpl_boundary((1,0,-1),s.Rational(1,2),region='0<x<1',
                              boundary='upper',max_terms=1)
    assert mp.mp.dps == original
