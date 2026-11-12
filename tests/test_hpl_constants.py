import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_constants import _remainder, evaluate_euler_sum, evaluate_mzv
from hypra.hpl_numeric import HPLNumericalConvergenceError
from hypra.symbolic import MZV, EulerSum


@pytest.mark.parametrize('indices,colors,expected', [
    ((1,), (-1,), -s.log(2)),
    ((2,), (-1,), -s.zeta(2)/2),
    ((1,1), (-1,-1), (s.log(2)**2-s.zeta(2))/2),
    ((2,2), (1,1), 3*s.zeta(4)/4),
    ((), (), s.S.One),
])
def test_exact_reducible_constants(indices, colors, expected):
    result = evaluate_euler_sum(indices, colors, dps=50)
    with mp.workdps(70):
        assert abs(result.value - mp.mpf(str(expected.evalf(70)))) < mp.mpf('1e-50')
    assert result.working_dps == 90


def test_mzv_unreduced_weight_six_against_independent_stuffle():
    # zeta(4,2)+zeta(2,4) = zeta(4)*zeta(2)-zeta(6).
    with mp.workdps(60):
        left = evaluate_mzv((4,2), dps=40).value + evaluate_mzv((2,4), dps=40).value
        assert abs(left - (mp.zeta(4)*mp.zeta(2)-mp.zeta(6))) < mp.mpf('1e-40')


def test_alternating_stuffle_with_distinct_indices():
    with mp.workdps(65):
        left = evaluate_euler_sum((1,2), (-1,-1), dps=45).value
        left += evaluate_euler_sum((2,1), (-1,-1), dps=45).value
        expected = mp.log(2)*mp.zeta(2)/2 - mp.zeta(3)
        assert abs(left - expected) < mp.mpf('1e-45')


@pytest.mark.parametrize('word', [(0,-1,1),(-1,0,1),(0,0,0,1,0,1)])
def test_symbolic_remainder_strictly_lowers_constant_weight(word):
    remainder = _remainder(word)
    assert all(sum(atom.indices) < len(word) for atom in remainder.atoms(EulerSum, MZV))


def test_reject_limits_divergence_and_invalid_inputs():
    with pytest.raises(UnsupportedExpressionError, match='weight'):
        evaluate_mzv((7,))
    with pytest.raises(ValueError):
        evaluate_euler_sum((1,), (1,))
    with pytest.raises(TypeError):
        evaluate_mzv('2,1')
    with pytest.raises(ValueError):
        evaluate_euler_sum((2.0,), (-1,))
    with pytest.raises(ValueError):
        evaluate_mzv((2,), dps=True)
    with pytest.raises(HPLNumericalConvergenceError):
        evaluate_euler_sum((2,1), (-1,-1), max_terms=1)


def test_global_precision_unchanged():
    before = mp.mp.dps
    evaluate_euler_sum((1,2), (-1,-1), dps=25)
    assert mp.mp.dps == before


def test_explicit_symbolic_constant_evalf_and_no_eager_constructor():
    mzv = MZV((4, 2))
    euler = EulerSum((1, 2), (-1, -1))
    assert isinstance(mzv, MZV)
    assert isinstance(euler, EulerSum)
    for atom, result in [
        (mzv, evaluate_mzv((4, 2), dps=55)),
        (euler, evaluate_euler_sum((1, 2), (-1, -1), dps=55)),
    ]:
        value = atom.evalf(50)
        assert isinstance(value, s.Float)
        with mp.workdps(70):
            assert abs(mp.mpf(str(value)) - result.value) < mp.mpf('1e-49')
    assert MZV((7,)).evalf(50) == MZV((7,))
    assert EulerSum((7,), (-1,)).evalf(50) == EulerSum((7,), (-1,))
