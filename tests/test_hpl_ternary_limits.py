import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_alternating import hpl_alternating_at_one
from hypra.hpl_limits import hpl_limit_at_one as limit
from hypra.symbolic import HPL

z=s.Symbol('z',positive=True)


def test_rational_poles_and_logarithmic_cancellations():
    shifted=HPL((-1,),z)-s.log(2)
    assert limit(shifted/(1-z),z) == -s.Rational(1,2)
    assert limit((shifted+(1-z)/2)/(1-z)**2,z) == -s.Rational(1,8)
    assert limit(HPL((1,),z)*(shifted/(1-z)+s.Rational(1,2)),z) == 0
    assert limit((HPL((-1,),z)**2-s.log(2)**2)/(1-z),z) == -s.log(2)


def test_convergent_ternary_endpoints():
    for word in ((-1,0),(-1,1),(0,-1,1),(-1,0,1)):
        assert s.simplify(limit(HPL(word,z),z)-hpl_alternating_at_one(word)) == 0


def test_divergent_and_unsupported_fail_closed():
    for expr in (HPL((-1,),z)/(1-z), HPL((1,-1),z),
                 1/HPL((-1,),z), HPL((-1,),z)/s.sin(z),
                 HPL((-1,),z)/(1-z)**17):
        with pytest.raises(UnsupportedExpressionError):
            limit(expr,z)
