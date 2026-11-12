"""Independent one-dimensional integrals for mixed-color endpoint constants."""

import mpmath as mp
import pytest

from hypra.hpl_constants import evaluate_euler_sum


@pytest.mark.parametrize("colors,sign,argument", [((-1, -1), -1, 1), ((-1, 1), -1, -1)])
def test_mixed_color_endpoint_against_independent_dilog_integral(colors, sign, argument):
    # H(-1,0,±1;1) has a single Li2 integrand. This reference uses neither
    # the Mobius transformation nor nested-sum endpoint recursion.
    result = evaluate_euler_sum((1, 2), colors, dps=25)
    context = mp.mp.clone()
    context.dps = 45
    reference = sign * context.quad(lambda t: context.polylog(2, argument * t) / (1 + t), [0, 1])
    assert abs(result.value - reference) < context.mpf("1e-24")
