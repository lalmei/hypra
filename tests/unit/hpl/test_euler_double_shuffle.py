"""Independent integral and direct-sum evidence for colored constant reductions."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl.alternating import reduce_euler_sum
from hypra.hpl.euler_reduction import euler_depth_two_reductions, euler_weight_reductions
from hypra.symbolic import MZV, EulerSum


@pytest.mark.parametrize("weight", [3, 5, 7])
def test_every_admissible_depth_two_sign_pattern_reduces(weight):
    table = euler_depth_two_reductions(weight)
    assert len(table) == 4 * (weight - 1) - 2
    assert all(not value.has(MZV, EulerSum, s.Float) for value in table.values())
    for indices, value in table.items():
        assert reduce_euler_sum(EulerSum(*zip(*indices, strict=True))) == value


def test_nielsen_s32_at_minus_one_from_double_shuffle():
    expected = s.pi**2 * s.zeta(3) / 12 - 29 * s.zeta(5) / 32
    assert reduce_euler_sum(EulerSum((4, 1), (-1, 1))) == expected
    # S32(-1)=1/4 integral_0^1 log(t)^2 log(1+t)^2/t dt follows
    # directly from Nielsen's defining logarithmic integral.
    context = mp.mp.clone()
    context.dps = 40
    reference = context.quad(
        lambda t: context.log(t) ** 2 * context.log1p(t) ** 2 / (4 * t), [0, 1]
    )
    assert abs(reference - context.mpf(str(expected.evalf(38)))) < context.mpf("1e-35")


@pytest.mark.parametrize("colors", list(product((-1, 1), repeat=2)))
@pytest.mark.parametrize("indices", [(4, 1), (4, 3), (5, 2)])
def test_color_mapping_against_independent_strict_sums(indices, colors):
    context = mp.mp.clone()
    context.dps = 30
    cutoff = 3000
    value, inner = context.zero, context.zero
    for n in range(1, cutoff + 1):
        value += colors[0] ** n * inner / n ** indices[0]
        inner += context.mpf(colors[1] ** n) / n ** indices[1]
    actual = context.mpf(str(reduce_euler_sum(EulerSum(indices, colors)).evalf(28)))
    # Absolute tail bound uses |inner|<=H_n<=1+log(n), regardless of colors.
    outer = indices[0]
    bound = (context.log(cutoff) + 1 + context.mpf(1) / (outer - 1)) / (
        (outer - 1) * cutoff ** (outer - 1)
    )
    assert abs(actual - value) < bound


@pytest.mark.parametrize("weight", range(2, 7))
def test_complete_colored_low_weight_space(weight):
    table = euler_weight_reductions(weight)
    assert len(table) == 4 * 3 ** (weight - 2)
    assert all(not value.has(MZV, s.Float) for value in table.values())
    remaining = set().union(*(value.atoms(EulerSum) for value in table.values()))
    assert remaining == ({EulerSum((5, 1), (-1, -1))} if weight == 6 else set())
    for indices, value in table.items():
        assert s.expand(reduce_euler_sum(EulerSum(*zip(*indices, strict=True))) - value) == 0


def test_even_weight_nielsen_constant_against_logarithmic_integral():
    result = reduce_euler_sum(EulerSum((3, 1), (-1, 1)))
    assert result.has(s.polylog(4, s.S.Half))
    assert not result.has(EulerSum)
    context = mp.mp.clone()
    context.dps = 40
    reference = context.quad(lambda t: -context.log(t) * context.log1p(t) ** 2 / (2 * t), [0, 1])
    assert abs(reference - context.mpf(str(result.evalf(38)))) < context.mpf("1e-35")


@pytest.mark.parametrize("colors", list(product((-1, 1), repeat=3)))
def test_depth_three_weight_five_against_strict_sum(colors):
    context = mp.mp.clone()
    context.dps = 30
    cutoff = 3000
    totals = [context.zero] * 3 + [context.one]
    for n in range(1, cutoff + 1):
        for depth, power in enumerate((3, 1, 1)):
            totals[depth] += colors[depth] ** n * totals[depth + 1] / n**power
    result = reduce_euler_sum(EulerSum((3, 1, 1), colors))
    actual = context.mpf(str(result.evalf(28)))
    bound = context.quad(lambda t: (1 + context.log(t)) ** 2 / (2 * t**3), [cutoff, context.inf])
    assert abs(actual - totals[0]) < bound


def test_higher_weight_primitive_is_preserved():
    value = EulerSum((6, 1, 1), (-1, 1, 1))
    assert reduce_euler_sum(value) == value


@pytest.mark.parametrize("solver", [euler_depth_two_reductions, euler_weight_reductions])
def test_cached_colored_reductions_cannot_be_mutated(solver):
    table = solver(5)
    word = ((4, -1), (1, 1))
    expected = table[word]
    with pytest.raises(TypeError):
        table[word] = s.S.Zero
    assert reduce_euler_sum(EulerSum((4, 1), (-1, 1))) == expected


@pytest.mark.parametrize("depth", [3, 7, 12, 24])
def test_all_unit_minus_minus_family_matches_half_polylog_integral(depth):
    result = reduce_euler_sum(EulerSum((1,) * depth, (-1, -1) + (1,) * (depth - 2)))
    if depth != 3:
        assert s.simplify(result + s.polylog(depth, s.S.Half)) == 0
    context = mp.mp.clone()
    context.dps = 40
    # u=e^-v removes endpoint logarithms and makes a rapidly decaying
    # integral that remains stable at the largest supported weight.
    reference = -context.quad(
        lambda v: v ** (depth - 1) / (2 * context.exp(v) - 1), [0, 1, depth, context.inf]
    ) / context.factorial(depth - 1)
    assert abs(reference - context.mpf(str(result.evalf(38)))) < context.mpf("1e-35")


@pytest.mark.parametrize("indices,sign", [((5, 1), 1), ((4, 1, 1), -1)])
def test_weight_six_reduction_against_independent_nielsen_integral(indices, sign):
    colors = (-1,) + (1,) * (len(indices) - 1)
    result = reduce_euler_sum(EulerSum(indices, colors))
    context = mp.mp.clone()
    context.dps = 38
    zeros, repeats = indices[0] - 1, len(indices)
    reference = (
        sign
        * context.quad(
            lambda t: (-context.log(t)) ** (zeros - 1) * context.log1p(t) ** repeats / t,
            [0, 1],
        )
        / (context.factorial(zeros - 1) * context.factorial(repeats))
    )
    assert abs(reference - context.mpf(str(result.evalf(32)))) < context.mpf("1e-28")
