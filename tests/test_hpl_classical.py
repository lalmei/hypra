from itertools import product

import pytest
import sympy as s

from hypra.hpl_classical import classical_hpl
from hypra.symbolic import HPL, reduce_known_hpl

z=s.Symbol('z',positive=True)

@pytest.mark.parametrize('word',list(product((-1,0,1),repeat=2)))
def test_all_weight_two_differential_and_origin(word):
    result=classical_hpl(HPL(word,z),z,domain='0<z<1')
    assert not result.has(HPL)
    difference=s.diff(result,z)-reduce_known_hpl(s.diff(HPL(word,z),z))
    difference=difference.xreplace({s.log((1-z)/2):s.log(1-z)-s.log(2),
                                    s.log((1+z)/2):s.log(1+z)-s.log(2)})
    assert s.simplify(s.expand_func(difference))==0
    if word!=(0,0):
        assert s.limit(result,z,0,dir='+')==0

def test_higher_mixed_words_remain_exact():
    value=HPL((1,-1,0,1),z)
    assert classical_hpl(value,z,domain='0<z<1')==value
