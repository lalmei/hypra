"""Independent low-weight and duality checks for the beta generating identity."""

import sympy as s

from hypra import MZV, reduce_mzv
from hypra.hpl.endpoints import _height_one, dual_indices


def test_known_height_one_values():
    assert s.simplify(reduce_mzv(MZV((3, 1))) - s.zeta(4) / 4) == 0
    assert s.simplify(reduce_mzv(MZV((4, 1))) - 2 * s.zeta(5) + s.zeta(2) * s.zeta(3)) == 0
    assert (
        s.simplify(reduce_mzv(MZV((5, 1))) - s.Rational(3, 4) * s.zeta(6) + s.zeta(3) ** 2 / 2) == 0
    )


def test_height_one_duality_and_depth_one():
    for first in range(2, 7):
        for depth in range(1, 5):
            assert s.simplify(_height_one(first, depth) - _height_one(depth + 1, first - 1)) == 0
            word = (first,) + (1,) * (depth - 1)
            assert reduce_mzv(MZV(word)) == reduce_mzv(MZV(dual_indices(word)))
        assert _height_one(first, 1) == s.zeta(first)


def test_unreduced_general_mzv_remains_exact():
    result = reduce_mzv(MZV((5, 3)))
    assert isinstance(result, MZV)


def test_euler_depth_two_and_repeated_indices_extend_to_weight_limit():
    first = 13
    expected = (
        first * s.zeta(first + 1) / 2
        - sum(s.zeta(first - k) * s.zeta(k + 1) for k in range(1, first - 1)) / 2
    )
    assert s.simplify(reduce_mzv(MZV((first, 1))) - expected) == 0
    for depth in (5, 8, 12):
        assert (
            s.simplify(
                reduce_mzv(MZV((2,) * depth)) - s.pi ** (2 * depth) / s.factorial(2 * depth + 1)
            )
            == 0
        )
    assert s.simplify(reduce_mzv(MZV((5, 5))) - (s.zeta(5) ** 2 - s.zeta(10)) / 2) == 0
    assert isinstance(reduce_mzv(MZV((2,) * 13)), MZV)


def test_zagier_family_against_generating_function_and_weight_eight_solver():
    from hypra.hpl.mzv_reduction import mzv_weight_reductions

    assert reduce_mzv(MZV((3, 1) * 2)) == mzv_weight_reductions(8)[(3, 1, 3, 1)]
    # Coefficients of (cosh(pi*t)-cos(pi*t))/(pi*t)^2 give the family.
    t = s.Symbol("t")
    generating = (s.cosh(s.pi * t) - s.cos(s.pi * t)) / (s.pi * t) ** 2
    series = generating.series(t, 0, 25).removeO().expand()
    for repeats in range(1, 7):
        assert reduce_mzv(MZV((3, 1) * repeats)) == series.coeff(t, 4 * repeats)


def test_prefixed_family_agrees_with_independent_double_shuffle_at_weight_six():
    from hypra.hpl.mzv_reduction import mzv_weight_reductions

    assert reduce_mzv(MZV((2, 1, 3))) == mzv_weight_reductions(6)[(2, 1, 3)]
    for repeats in range(1, 6):
        indices = (2,) + (1, 3) * repeats
        result = reduce_mzv(MZV(indices))
        assert not result.has(MZV)
        assert result.evalf(30).is_positive
        assert result == reduce_mzv(MZV(dual_indices(indices)))
