from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl_classical import classical_hpl
from hypra.hpl_numeric import evaluate_hpl
from hypra.symbolic import HPL

Z = s.Symbol("z", positive=True)
DOMAIN = {"domain": "0<z<1"}


@pytest.mark.parametrize("word", list(product((-1, 1), repeat=3)))
def test_mixed_two_letter_cubic_values(word):
    result = classical_hpl(HPL(word, Z), Z, **DOMAIN)
    assert not result.has(HPL, s.I)
    context = mp.mp.clone()
    context.dps = 40
    for point in (s.Rational(1, 3), s.Rational(3, 4)):
        actual = context.mpf(str(s.N(result.subs(Z, point), 38)))
        reference = evaluate_hpl(word, point, dps=30).value
        assert abs(actual - reference) < context.mpf("1e-29")


@pytest.mark.parametrize("word", [(-1, 1, 1), (1, -1, -1)])
def test_integral_primitive_derivative_and_origin(word):
    result = classical_hpl(HPL(word, Z), Z, **DOMAIN)
    tail = classical_hpl(HPL(word[1:], Z), Z, **DOMAIN)
    kernel = 1 / (1 - Z) if word[0] == 1 else 1 / (1 + Z)
    # Every log argument is positive on the declared interval.
    residual = s.expand_func(s.diff(result, Z) - kernel * tail)
    u = (1 - Z) / 2
    # Euler reflection is independently established by derivative and value at0.
    residual = residual.subs(
        s.polylog(2, (1 + Z) / 2), s.pi**2 / 6 - s.polylog(2, u) - s.log(u) * s.log(1 - u)
    )
    residual = s.expand_log(residual, force=True)
    assert s.simplify(residual) == 0
    assert s.limit(result, Z, 0, dir="+") == 0


def test_true_ternary_primitives_preserved():
    for word in [(-1, 0, 1), (-1, 1, 0)]:
        atom = HPL(word, Z)
        assert classical_hpl(atom, Z, **DOMAIN) == atom
