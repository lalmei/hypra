import pytest
import sympy as s

from hypra import HPL, UnsupportedExpressionError, integrate_from_zero, primitive_hpl
from hypra.hpl_algebra import shuffle_expand

z=s.Symbol('z',positive=True)

@pytest.mark.parametrize('expression',[
    HPL((0,),z)/z, HPL((1,),z)/(1-z),
    z**2*HPL((0,1),z)/(1-z**2)**2,
    HPL((1,-1),z)/(1+z)**3,
    HPL((1,),z)*HPL((-1,),z)/z,
    s.log(z)/(z**2),
])
def test_exact_derivative(expression):
    answer=primitive_hpl(expression,z)
    assert s.simplify(shuffle_expand(s.diff(answer,z)-expression,z))==0

def test_divergent_definite_not_conflated_with_primitive():
    assert primitive_hpl(1/z,z)==HPL((0,),z)
    with pytest.raises(UnsupportedExpressionError):
        integrate_from_zero(1/z,z)
