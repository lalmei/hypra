"""Proved identities and independent finite-sum checks for MZV reduction."""

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl.endpoints import reduce_mzv
from hypra.hpl.mzv_reduction import mzv_weight_reductions
from hypra.symbolic import MZV


@pytest.mark.parametrize("weight", range(2, 9))
def test_complete_admissible_weight_space_and_chosen_primitive(weight):
    reductions = mzv_weight_reductions(weight)
    assert len(reductions) == 2 ** (weight - 2)
    remaining = set().union(*(value.atoms(MZV) for value in reductions.values()))
    assert remaining == ({MZV((5, 3))} if weight == 8 else set())
    assert all(reduce_mzv(value) == value for value in reductions.values())


@pytest.mark.parametrize("depth", range(1, 5))
def test_all_twos_from_euler_sine_product(depth):
    # Comparing coefficients in sinh(pi*sqrt(t))/(pi*sqrt(t)) proves this
    # independently of the double-shuffle equations used by the reducer.
    expected = s.pi ** (2 * depth) / s.factorial(2 * depth + 1)
    assert s.simplify(reduce_mzv(MZV((2,) * depth)) - expected) == 0


def test_double_zeta_stuffle_and_known_euler_reduction():
    assert reduce_mzv(MZV((3, 2))) == s.pi**2 * s.zeta(3) / 2 - 11 * s.zeta(5) / 2
    assert s.simplify(reduce_mzv(MZV((3, 3))) - (s.zeta(3) ** 2 - s.zeta(6)) / 2) == 0
    assert s.simplify(reduce_mzv(MZV((4, 2))) - s.zeta(3) ** 2 + 4 * s.zeta(6) / 3) == 0


@pytest.mark.parametrize("indices", [(3, 2), (4, 2), (5, 2), (4, 3), (4, 2, 1)])
def test_reduced_values_against_independent_strict_sums(indices):
    context = mp.mp.clone()
    context.dps = 30
    cutoff = 2500
    tails = [context.zero] * len(indices) + [context.one]
    for n in range(1, cutoff + 1):
        # Ascending depth means each inner tail still has the n-1 cutoff.
        for depth, power in enumerate(indices):
            tails[depth] += tails[depth + 1] / n**power
    actual = context.mpf(str(reduce_mzv(MZV(indices)).evalf(28)))
    # Inner sums are bounded by H_n^(depth-1)/(depth-1)!. An integral of
    # the positive envelope gives a rigorous bound on the uncomputed tail.
    outer, depth = indices[0], len(indices)
    bound = context.quad(
        lambda x: (1 + context.log(x)) ** (depth - 1) / (context.factorial(depth - 1) * x**outer),
        [cutoff, context.inf],
    )
    assert 0 < actual - tails[0] < bound


def test_higher_weights_remain_exact():
    assert isinstance(reduce_mzv(MZV((7, 3))), MZV)
    with pytest.raises(ValueError):
        mzv_weight_reductions(9)


def test_cached_mzv_reductions_cannot_be_mutated():
    table = mzv_weight_reductions(5)
    expected = table[(3, 2)]
    with pytest.raises(TypeError):
        table[(3, 2)] = s.S.Zero
    assert reduce_mzv(MZV((3, 2))) == expected
