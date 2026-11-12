import pytest
import sympy as s

from hypra import HPL, UnsupportedExpressionError
from hypra.hpl_arguments import transform_hpl_arguments
from hypra.symbolic import reduce_known_hpl

z=s.Symbol('z',positive=True)

def test_mixed_expression_unification():
    original=HPL((1,),1-z)+2*HPL((0,),1/z)+HPL((-1,),(1-z)/(1+z))
    result=transform_hpl_arguments(original,z,domain='0<z<1',boundary='lower')
    assert all(h.argument==z for h in result.atoms(HPL))
    assert s.simplify(reduce_known_hpl(result)-(-3*s.log(z)+s.log(2)-s.log(1+z)))==0

@pytest.mark.parametrize('argument',[-z,1/z,z/(z-1)])
def test_cut_selection_required(argument):
    with pytest.raises(UnsupportedExpressionError,match='boundary'):
        transform_hpl_arguments(HPL((0,),argument),z,domain='0<z<1')

def test_unknown_variable_argument_rejected_and_constants_preserved():
    with pytest.raises(UnsupportedExpressionError,match='Unsupported HPL argument'):
        transform_hpl_arguments(HPL((1,),s.sqrt(z)),z,domain='0<z<1')
    constant=HPL((0,1),s.I)
    assert transform_hpl_arguments(constant,z,domain='0<z<1')==constant
