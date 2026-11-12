from itertools import product

import pytest
import sympy as s

from hypra.hpl_notation import (
    compress_indices,
    expand_indices,
    hpl_pm,
    letters_to_pm,
    pm_to_letters,
)
from hypra.symbolic import HPL


def test_all_short_word_roundtrips():
    for n in range(5):
        for word in product((-1,0,1),repeat=n):
            assert expand_indices(compress_indices(word))==word
            result={}
            for pm,coefficient in letters_to_pm(word).items():
                for letters,c in pm_to_letters(pm).items():
                    result[letters]=result.get(letters,0)+coefficient*c
            assert {w:c for w,c in result.items() if c}=={word:1}

def test_kernel_conventions():
    z=s.Symbol('z')
    assert hpl_pm(('plus',),z)==HPL((1,),z)+HPL((-1,),z)
    assert hpl_pm(('minus',),z)==HPL((1,),z)-HPL((-1,),z)
    assert s.simplify(s.diff(hpl_pm(('plus',),z),z)-2/(1-z**2))==0
    assert s.simplify(s.diff(hpl_pm(('minus',),z),z)-2*z/(1-z**2))==0
    assert expand_indices((2,-3,0))==(0,1,0,0,-1,0)

@pytest.mark.parametrize('word',[(True,),(0.0,),('bad',)])
def test_bad_pm_letters(word):
    with pytest.raises((ValueError,TypeError)):
        pm_to_letters(word)
