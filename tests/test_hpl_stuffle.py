import itertools

import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_stuffle import stuffle_expand_constants as expand
from hypra.symbolic import MZV, EulerSum


def cutoff(atom, nmax):
    if not isinstance(atom, (MZV, EulerSum)):
        return atom
    colors = (1,)*len(atom.indices) if isinstance(atom, MZV) else atom.colors
    return sum(s.prod(s.Integer(c)**n/s.Integer(n)**k for c,n,k in zip(colors, ns, atom.indices))
               for ns in itertools.combinations(range(nmax, 0, -1), len(atom.indices)))


def test_mzv_stuffle_and_repeated_multiplicity():
    assert expand(MZV((2,))**2) == 2*MZV((2,2))+MZV((4,))
    assert expand(MZV((2,))*MZV((3,))) == MZV((2,3))+MZV((3,2))+MZV((5,))


@pytest.mark.parametrize('left,right', [
    (EulerSum((1,),(-1,)), EulerSum((1,),(-1,))),
    (MZV((2,1)), EulerSum((1,),(-1,))),
    (EulerSum((1,2),(-1,-1)), EulerSum((2,),(1,))),
])
def test_exact_finite_cutoff_identity(left,right):
    result = expand(left*right)
    evaluated = result.xreplace({a:cutoff(a,8) for a in result.atoms(MZV,EulerSum)})
    assert evaluated == cutoff(left,8)*cutoff(right,8)


def test_polynomial_and_fail_closed():
    a=MZV((2,));x=s.Symbol('x')
    assert s.expand(expand((x+a)**2)-x*x-2*x*a-2*MZV((2,2))-MZV((4,))) == 0
    for expr in (1/a,s.sin(a),a**9,MZV((7,))**2,s.Float(0.1)*a):
        with pytest.raises(UnsupportedExpressionError):
            expand(expr)
