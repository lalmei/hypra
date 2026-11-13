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


@pytest.mark.parametrize("name,word", [
    ("HPLs6", (-5, -1)),
    ("HPLs7a", (-5, 1, 1)),
    ("HPLs7b", (5, -1, -1)),
    ("HPLs8a", (5, 3)),
    ("HPLs8b", (-7, -1)),
    ("HPLs8c", (-5, -1, -1, -1)),
    ("HPLs8d", (-5, -1, 1, 1)),
])
def test_named_constant_against_independent_weak_sum(name, word):
    import hypra

    context = mp.mp.clone()
    context.dps = 25
    sums = [context.zero] * len(word) + [context.one]
    for index in range(1, 2001):
        for j in reversed(range(len(word))):
            power = word[j]
            sign = -1 if power < 0 and index % 2 else 1
            sums[j] += sign * sums[j+1] / index**abs(power)
    value = context.mpf(str(getattr(hypra, name).evalf(16)))
    # For n > 2000, |inner sums| <= (1+log n)**3. Integrating
    # that envelope against n**-5 bounds the omitted terms by 1.1e-11.
    assert abs(value-sums[0]) < context.mpf("1.1e-11")
