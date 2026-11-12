import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_continue import continue_hpl
from hypra.hpl_numeric import HPLNumericalConvergenceError


def call(word,z,**kwargs):
    return continue_hpl(word,z,domain='principal-cut-plane',path='radial',**kwargs)


@pytest.mark.parametrize('word,kind', [((1,),'log'),((0,1),'li2'),((0,0,1),'li3'),
                                      ((0,),'logz'),((-1,),'logplus')])
def test_imaginary_beyond_unit_disk(word,kind):
    result=call(word,3*s.I,dps=25)
    with mp.workdps(60):
        z=3*mp.j
        expected={'log':-mp.log(1-z),'li2':mp.polylog(2,z),'li3':mp.polylog(3,z),
                  'logz':mp.log(z),'logplus':mp.log(1+z)}[kind]
        assert abs(result.value-expected)<mp.mpf('1e-25')
        assert result.tail_bound<mp.mpf('1e-25')
        assert result.stability_gap<mp.mpf('1e-25')
    assert result.steps>1


def test_suffix_trailing_log_and_shuffle_identity():
    a,b=call((1,0),2*s.I,dps=20),call((0,1),2*s.I,dps=20)
    with mp.workdps(50):
        assert abs(a.value+b.value+mp.log(1-2*mp.j)*mp.log(2*mp.j))<mp.mpf('1e-20')


@pytest.mark.parametrize('z',[0,1,-1,2,-2])
def test_cuts_endpoints_rejected(z):
    with pytest.raises(UnsupportedExpressionError):
        call((1,),z,dps=15)


def test_domain_and_resource_requirements():
    with pytest.raises(TypeError):
        continue_hpl((1,),2*s.I)
    with pytest.raises(UnsupportedExpressionError):
        continue_hpl((1,),2*s.I,domain='automatic',path='radial')
    with pytest.raises(HPLNumericalConvergenceError):
        call((0,1),3*s.I,dps=20,max_steps=1)
    with pytest.raises(HPLNumericalConvergenceError):
        call((1,),3*s.I,dps=20,max_degree=2)


def test_lower_half_plane_and_initial_disk_shortcut():
    lower=call((1,),-2*s.I,dps=20)
    with mp.workdps(50):
        assert abs(lower.value+mp.log(1+2*mp.j))<mp.mpf('1e-20')
    inside=call((0,1),s.Rational(1,8),dps=20)
    assert inside.steps==0
    assert call((),2*s.I,dps=15).value==1
