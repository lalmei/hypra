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
