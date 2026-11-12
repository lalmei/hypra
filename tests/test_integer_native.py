"""Actual native pF(p-1) coefficients vs independent finite Pochhammer quotients.

The reference is generated directly for each argument power; it contains no
stored hypergeometric expansion coefficients and invokes no Wolfram backend.
"""

import pytest
import sympy as s

from hypra import NativeConfig, expand_epsilon
from hypra.hpl_algebra import hpl_series

E, Z = s.symbols("epsilon z")
DOMAIN = s.And(Z > 0, Z < 1)
NO_CACHE = NativeConfig(use_cache=False)


def rational_laurent(expression, through):
    """Exact formal division of a rational regulator expression at zero."""
    expression = s.cancel(expression)
    if expression == 0:
        return {}
    numerator, denominator = s.fraction(expression)
    numerator, denominator = s.Poly(numerator, E), s.Poly(denominator, E)
    numerator_valuation = min(power[0] for power, _ in numerator.terms())
    denominator_valuation = min(power[0] for power, _ in denominator.terms())
    valuation = numerator_valuation - denominator_valuation
    regular_denominator = [
        denominator.nth(denominator_valuation + j) for j in range(through - valuation + 1)
    ]
    if through < valuation:
        return {}
    result = []
    for k in range(through - valuation + 1):
        previous = sum(regular_denominator[j] * result[k - j] for j in range(1, k + 1))
        result.append(
            s.cancel((numerator.nth(numerator_valuation + k) - previous) / regular_denominator[0])
        )
    return {valuation + k: value for k, value in enumerate(result)}


def direct_argument_coefficients(upper, lower, through, argument_order):
    """Finite z-series reference from the defining generalized hypergeometric series."""
    result = {}
    for n in range(argument_order):
        coefficient = s.prod(s.rf(a, n) for a in upper) / (
            s.factorial(n) * s.prod(s.rf(b, n) for b in lower)
        )
        for power, value in rational_laurent(coefficient, through).items():
            result[power] = result.get(power, s.S.Zero) + value * Z**n
    return result


CASES = [
    pytest.param([2 + E, 1 - E], [3 + 2 * E], 2, 4, id="2f1-positive-shifts"),
    pytest.param(
        [1 + E, 2 - 2 * E, 1 + 2 * E], [2 + 3 * E, 1 - 3 * E], 2, 4, id="3f2-positive-shifts"
    ),
    pytest.param(
        [1 + E, 1 - E, 1 + 2 * E, 1 - 2 * E],
        [1 + 3 * E, 1 - 3 * E, 1 + 4 * E],
        2,
        4,
        id="4f3-canonical",
    ),
    pytest.param([E, 1 + 2 * E], [2 - E], 3, 5, id="2f1-zero-upper"),
    pytest.param([E, 2 * E, 1 - E], [1 + E, 2 - 3 * E], 3, 4, id="3f2-two-zero-uppers"),
    pytest.param([-1 + E, 1 + 2 * E], [1 - 2 * E], 2, 5, id="2f1-negative-moving-upper"),
    pytest.param(
        [-1 + E, 1 + E, 1 - E], [1 - 2 * E, 1 + 2 * E], 2, 4, id="3f2-negative-moving-upper"
    ),
    pytest.param([1 + E, 1 - 2 * E], [E], 2, 4, id="2f1-zero-moving-lower"),
    pytest.param([1 + E, 2 - E], [-1 + 2 * E], 2, 5, id="2f1-negative-moving-lower"),
    pytest.param([1, 1 - E, 1 + E], [E, 2 * E], 0, 3, id="3f2-two-singular-lowers"),
    pytest.param([-2, 1 + E], [2 - E], 3, 4, id="2f1-fixed-terminating-upper"),
    pytest.param([1, 1 + E, 1 - E], [E, 2 + E], 1, 4, id="3f2-singular-and-positive-lower"),
]


@pytest.mark.parametrize("upper,lower,through,argument_order", CASES)
def test_integer_native_against_direct_pochhammer_coefficients(
    upper, lower, through, argument_order
):
    expression = s.hyper(upper, lower, Z, evaluate=False)
    actual = expand_epsilon(expression, E, through, DOMAIN, config=NO_CACHE)
    reference = direct_argument_coefficients(upper, lower, through, argument_order)
    lower_bound = min(reference)
    for power in range(lower_bound, through + 1):
        calculated = hpl_series(actual.coefficient(power), Z, argument_order)
        residual = s.simplify(s.expand(calculated - reference.get(power, s.S.Zero)))
        assert residual == 0, (expression, power, residual)
    # No additional spurious Laurent poles may be introduced by the shifts.
    for power in range(-len(lower) - 1, lower_bound):
        assert actual.coefficient(power) == 0
    assert actual.provenance["backend"]["backend"] == "native-sympy"
    assert actual.provenance["injected_backend"] is False


MATCHED_HALF_CASES = [
    pytest.param(
        [1 + E, s.Rational(1, 2) - E], [s.Rational(1, 2) + 2 * E], 2, 4, id="half-canonical-2f1"
    ),
    pytest.param(
        [2 + E, s.Rational(3, 2) - E], [s.Rational(1, 2) + E], 1, 4, id="half-upper-shifts"
    ),
    pytest.param(
        [1 + E, s.Rational(1, 2) - E], [s.Rational(3, 2) + E], 2, 4, id="half-lower-raising"
    ),
    pytest.param(
        [1 + E, s.Rational(1, 2) - E, s.Rational(1, 2) + 2 * E],
        [s.Rational(1, 2) + 3 * E, s.Rational(1, 2) - 2 * E],
        1,
        4,
        id="two-half-pairs-3f2",
    ),
    pytest.param(
        [E, s.Rational(1, 2) + E, 1 + 2 * E],
        [s.Rational(1, 2) - E, 2 + E],
        2,
        4,
        id="half-zero-upper-3f2",
    ),
    pytest.param(
        [-s.Rational(1, 2) + E, 1 + E],
        [s.Rational(1, 2) + 2 * E],
        1,
        4,
        id="half-negative-moving-upper",
    ),
    pytest.param(
        [1, s.Rational(1, 2) + E], [-s.Rational(1, 2) + 2 * E], 2, 4, id="half-negative-lower"
    ),
]


@pytest.mark.parametrize("upper,lower,through,argument_order", MATCHED_HALF_CASES)
def test_matched_half_native_against_direct_pochhammer_coefficients(
    upper, lower, through, argument_order
):
    expression = s.hyper(upper, lower, Z, evaluate=False)
    actual = expand_epsilon(expression, E, through, DOMAIN, config=NO_CACHE)
    reference = direct_argument_coefficients(upper, lower, through, argument_order)
    t = s.Symbol("t", positive=True)
    for power in range(min(reference), through + 1):
        coefficient = s.expand_log(actual.coefficient(power).subs(Z, t**2), force=False)
        calculated = hpl_series(coefficient, t, 2 * argument_order)
        residual = s.simplify(s.expand(calculated - reference.get(power, s.S.Zero).subs(Z, t**2)))
        assert residual == 0, (expression, power, residual)
    assert actual.provenance["backend"]["backend"] == "native-sympy"


@pytest.mark.parametrize(
    "upper,lower",
    [
        pytest.param([1 + E, E], [s.Rational(1, 2) + 2 * E], id="shifted-pattern-A"),
        pytest.param([s.Rational(1, 2) + E, 1 + E], [1 + 2 * E], id="shifted-pattern-B"),
        pytest.param(
            [s.Rational(3, 2) + E, s.Rational(1, 2) - E],
            [s.Rational(1, 2) + E],
            id="shifted-pattern-D",
        ),
    ],
)
def test_shifted_twof1_rationalized_systems_against_pochhammer(upper, lower):
    from hypra.symbolic import HPL

    expression = s.hyper(upper, lower, Z, evaluate=False)
    actual = expand_epsilon(expression, E, 1, DOMAIN, config=NO_CACHE)
    reference = direct_argument_coefficients(upper, lower, 1, 4)
    t = s.Symbol("local_t", positive=True)
    for power in range(2):
        coefficient = actual.coefficient(power)
        # All these rationalizing arguments vanish at z=0. Retain more local
        # HPL terms than needed for z^3, including their rational prefactors.
        # Complex substitution uses the actual principal branch from the result.
        coefficient = coefficient.xreplace(
            {
                atom: hpl_series(HPL(atom.word, t), t, 12).subs(t, atom.argument)
                for atom in coefficient.atoms(HPL)
            }
        )
        calculated = s.series(coefficient, Z, 0, 4).removeO()
        assert s.simplify(s.expand(calculated - reference[power])) == 0
