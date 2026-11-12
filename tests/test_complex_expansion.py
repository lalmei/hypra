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
