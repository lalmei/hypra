import pytest
import sympy as s

from hypra.hpl_algebra import hpl_series
from hypra.hpl_reduce import extract_hpl_logs
from hypra.symbolic import HPL

z=s.Symbol('z',positive=True)

@pytest.mark.parametrize('word',[(1,0),(1,-1,0),(1,1,-1,0,0),(0,0),(-1,-1),(1,1,1),(0,1,0),(1,0,1)])
def test_shuffle_log_extraction(word):
    original=HPL(word,z)
    result=extract_hpl_logs(original)
    assert s.simplify(hpl_series(result-original,z,6)) == 0
    for atom in result.atoms(HPL):
        assert atom.word[0] != 1 and atom.word[-1] != 0

def test_expected_weight_two():
    assert s.simplify(extract_hpl_logs(HPL((1,0),z)) + s.log(z)*s.log(1-z)+s.polylog(2,z)) == 0

def test_flags_and_validation():
    assert extract_hpl_logs(HPL((1,0),z),leading_ones=False,trailing_zeros=False)==HPL((1,0),z)
    with pytest.raises(TypeError):
        extract_hpl_logs('text')


def test_lyndon_factorization_and_basis():
    from itertools import product

    from hypra.hpl_reduce import lyndon_factors, lyndon_reduce
    for weight in range(1,4):
        for word in product((-1,0,1),repeat=weight):
            factors=lyndon_factors(word)
            assert sum(factors,())==word
            assert list(factors)==sorted(factors,reverse=True)
            result=lyndon_reduce(HPL(word,z))
            assert s.simplify(hpl_series(result-HPL(word,z),z,5))==0
            for atom in result.atoms(HPL):
                assert len(lyndon_factors(atom.word))==1


def test_lyndon_repeated_and_limits():
    from hypra.errors import UnsupportedExpressionError
    from hypra.hpl_reduce import lyndon_reduce
    assert lyndon_reduce(HPL((1,1,1),z))==HPL((1,),z)**3/6
    with pytest.raises(UnsupportedExpressionError,match='eight'):
        lyndon_reduce(HPL((0,)*8+(1,),z))
