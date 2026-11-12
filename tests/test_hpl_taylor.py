import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_taylor import hpl_taylor
from hypra.symbolic import HPL

z=s.Symbol('z');h=s.Symbol('h');c=s.Rational(1,3)


@pytest.mark.parametrize('word', [(),(0,),(1,),(-1,),(0,1),(-1,0,1),(1,0,-1)])
def test_taylor_coefficients_equal_exact_derivatives(word):
    result=hpl_taylor(word,c,z,through=4).removeO().subs(z,c+h).expand()
    function=HPL(word,z) if word else s.S.One
    for n in range(5):
        assert s.simplify(result.coeff(h,n)-function.subs(z,c)/s.factorial(n)) == 0
        function=s.diff(function,z)


def test_sympy_series_at_shifted_point():
    actual=HPL((0,1),c+h).series(h,0,4).removeO()
    expected=hpl_taylor((0,1),c,z,through=3).removeO().subs(z,c+h)
    assert s.expand(actual-expected) == 0
    actual=HPL((-1,0,1),c+h*h).series(h,0,5).removeO()
    expected=hpl_taylor((-1,0,1),c,z,through=2).removeO().subs(z,c+h*h)
    assert s.expand(actual-expected) == 0


@pytest.mark.parametrize('center', [0,1,-1,s.Float(.5),s.I,s.Symbol('c')])
def test_invalid_center_fails(center):
    with pytest.raises(UnsupportedExpressionError):
        hpl_taylor((0,1),center,z)
