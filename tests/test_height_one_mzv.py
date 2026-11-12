"""Independent low-weight and duality checks for the beta generating identity."""
import sympy as s

from hypra import MZV, reduce_mzv
from hypra.hpl_endpoints import _height_one, dual_indices


def test_known_height_one_values():
    assert s.simplify(reduce_mzv(MZV((3, 1))) - s.zeta(4)/4) == 0
    assert s.simplify(reduce_mzv(MZV((4, 1))) - 2*s.zeta(5) + s.zeta(2)*s.zeta(3)) == 0
    assert s.simplify(reduce_mzv(MZV((5, 1))) - s.Rational(3,4)*s.zeta(6) + s.zeta(3)**2/2) == 0


def test_height_one_duality_and_depth_one():
    for first in range(2, 7):
        for depth in range(1, 5):
            assert s.simplify(_height_one(first, depth) - _height_one(depth+1, first-1)) == 0
            word = (first,) + (1,)*(depth-1)
            assert reduce_mzv(MZV(word)) == reduce_mzv(MZV(dual_indices(word)))
        assert _height_one(first, 1) == s.zeta(first)


def test_unreduced_general_mzv_remains_exact():
    result = reduce_mzv(MZV((5, 3)))
    assert isinstance(result, MZV)
