"""Generic contiguous operators against independently generated coefficients."""

import pytest
import sympy as s

from hypra.contiguous import expand_contiguous, plan_contiguous
from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import hpl_series
from hypra.nested_sums import canonical_hyper_series
from hypra.planning import laurent_coefficient
from hypra.symbolic import reduce_known_hpl

E, Z = s.symbols("epsilon z")


def canonical_builder(bases, slopes, z, epsilon, depth):
    upper_bases, lower_bases = bases
    upper_slopes, lower_slopes = slopes
    assert all(base == 1 for base in lower_bases)
    return canonical_hyper_series(
        list(zip(upper_bases, upper_slopes, strict=True)), lower_slopes, z, epsilon, depth
    )


def _direct_reference(upper, lower, through, argument_order):
    result = s.S.Zero
    for n in range(argument_order):
        coefficient = s.prod(s.rf(base + slope * E, n) for base, slope in upper) / (
            s.factorial(n) * s.prod(s.rf(base + slope * E, n) for base, slope in lower)
        )
        result += Z**n * s.series(s.cancel(coefficient), E, 0, through + 1).removeO()
    return s.expand(result)


@pytest.mark.parametrize(
    "upper,lower,canonical_upper,through",
    [
        pytest.param([(-1, 1), (1, 2)], [(1, -2)], [0, 1], 2, id="negative-upper-2f1"),
        pytest.param(
            [(-2, 1), (1, 1), (2, -1)],
            [(1, -2), (2, 1)],
            [0, 1, 1],
            1,
            id="negative-upper-3f2-degenerate-shift",
        ),
        pytest.param([(1, 1), (1, -2)], [(0, 1)], [1, 1], 1, id="singular-lower"),
        pytest.param([(1, 1), (1, -2)], [(2, 3)], [1, 1], 1, id="positive-lower"),
    ],
)
def test_contiguous_shifts_against_direct_pochhammer(upper, lower, canonical_upper, through):
    result = expand_contiguous(
        upper, lower, canonical_upper, [1] * len(lower), Z, E, through, canonical_builder
    )
    reference = _direct_reference(upper, lower, through, 4)
    pole_bound = sum(1 for base, _ in lower if base <= 0)
    for power in range(-pole_bound, through + 1):
        actual = hpl_series(laurent_coefficient(result, E, power), Z, 4)
        expected = laurent_coefficient(reference, E, power)
        assert s.simplify(s.expand(actual - expected)) == 0


def test_one_delta_removable_limit_recovers_fixed_parameters():
    upper, lower = [(2, 0), (1, 0)], [(3, 0)]
    plan = plan_contiguous(upper, lower, [1, 1], [1])
    assert plan.perturbation_used
    assert plan.extra_order == 2
    result = expand_contiguous(upper, lower, [1, 1], [1], Z, E, 0, canonical_builder)
    expected = -2 * (Z + s.log(1 - Z)) / Z**2
    assert s.simplify(reduce_known_hpl(result) - expected) == 0
    assert result.free_symbols == {Z}


def test_nonzero_divisor_does_not_introduce_slope_perturbation():
    plan = plan_contiguous([(-1, 1), (1, 2)], [(1, -2)], [1, 1], [1])
    assert plan.extra_order == 1
    assert plan.perturbation_used is False


def test_resource_limits_and_fixed_lower_pole_reject():
    with pytest.raises(UnsupportedExpressionError, match="shift count"):
        plan_contiguous([(100, 1), (1, 2)], [(1, -2)], [1, 1], [1])
    with pytest.raises(UnsupportedExpressionError, match="Fixed singular"):
        plan_contiguous([(1, 1), (1, 2)], [(0, 0)], [1, 1], [1])
    with pytest.raises(UnsupportedExpressionError, match="depth"):
        expand_contiguous(
            [(1, 1), (1, -2)],
            [(2, 3)],
            [1, 1],
            [1],
            Z,
            E,
            2,
            canonical_builder,
            max_internal_order=3,
        )
