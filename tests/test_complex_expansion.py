import pytest
import sympy as s

from hypra import NativeConfig, UnsupportedExpressionError, expand_epsilon

e=s.Symbol('eps')

@pytest.mark.parametrize('z',[s.I/3,(1+s.I)/4,(1-s.I)/4])
def test_explicit_complex_disk_exact_p_family(z):
    result=expand_epsilon(s.hyper([1,-e],[1-2*e],z),e,2,
                          config=NativeConfig(argument_domain='principal-disk',use_cache=False))
    expected=1+e*s.log(1-z)+e**2*(-2*s.polylog(2,z)-s.log(1-z)**2/2)
    assert s.simplify(result.polynomial-expected)==0
    assert result.provenance['backend']['argument_domain']=='principal-disk'

@pytest.mark.parametrize('z',[-s.Rational(1,2),2*s.I,s.Symbol('z')])
def test_disk_contract_rejects_unproved_or_outside_points(z):
    with pytest.raises(UnsupportedExpressionError,match='Principal-disk'):
        expand_epsilon(s.hyper([1,-e],[1-2*e],z),e,1,
                       config=NativeConfig(argument_domain='principal-disk'))

def test_default_does_not_widen_complex_domain():
    with pytest.raises(UnsupportedExpressionError,match='0 <'):
        expand_epsilon(s.hyper([1,-e],[1-2*e],s.I/3),e,1)


@pytest.mark.parametrize('z',[3*s.I,2+s.I,2-s.I])
def test_explicit_cut_plane_closed_forms(z):
    result=expand_epsilon(s.hyper([1,-e],[1-2*e],z),e,2,
                          config=NativeConfig(argument_domain='principal-cut-plane',use_cache=False))
    expected=1+e*s.log(1-z)+e**2*(-2*s.polylog(2,z)-s.log(1-z)**2/2)
    assert s.simplify(result.polynomial-expected)==0
    # Secondary numerical check at a small nonzero regulator, not the derivation.
    numeric=s.N(result.polynomial.subs(e,s.Rational(1,100000)),35)
    oracle=s.N(s.hyper([1,-s.Rational(1,100000)],[1-s.Rational(2,100000)],z),35)
    assert abs(complex(numeric-oracle))<1e-12


@pytest.mark.parametrize('z',[-1,s.Rational(-1,2),1,2,s.Symbol('z')])
def test_cut_plane_requires_proved_nonboundary_point(z):
    with pytest.raises(UnsupportedExpressionError):
        expand_epsilon(s.hyper([1,-e],[1-2*e],z),e,1,
                       config=NativeConfig(argument_domain='principal-cut-plane'))
