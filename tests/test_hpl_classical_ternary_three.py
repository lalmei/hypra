from itertools import permutations

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl_alternating import hpl_alternating_at_one
from hypra.hpl_classical import classical_hpl
from hypra.hpl_mobius import mobius_hpl
from hypra.hpl_numeric import evaluate_hpl
from hypra.hpl_reduce import lyndon_reduce
from hypra.symbolic import HPL

Z = s.Symbol("z", positive=True)
DOMAIN = {"domain": "0<z<1"}


@pytest.mark.parametrize("word", list(permutations((-1, 0, 1))))
def test_every_genuinely_ternary_cubic_numerically(word):
    result = classical_hpl(HPL(word, Z), Z, **DOMAIN)
    assert not result.has(HPL, s.I)
    context = mp.mp.clone()
    context.dps = 40
    for point in (s.Rational(1, 3), s.Rational(3, 4)):
        numerical = context.mpf(str(s.N(result.subs(Z, point), 38)))
        reference = evaluate_hpl(word, point, dps=30).value
        assert abs(numerical - reference) < context.mpf("1e-29")


def test_first_primitive_exact_differential_and_origin():
    result = classical_hpl(HPL((-1, 0, 1), Z), Z, **DOMAIN)
    # The constituent real primitives have separately checked derivatives.
    # Verify the square-log construction before classical substitution here.
    primitive = (
        HPL((0, 1, 1), Z**2) / 2
        - HPL((0, -1, -1), Z)
        - HPL((0, 1, 1), Z)
        + s.log(1 + Z) * s.polylog(2, Z)
    )
    from hypra.symbolic import reduce_known_hpl

    derivative = reduce_known_hpl(s.diff(primitive, Z) - s.polylog(2, Z) / (1 + Z))
    derivative = s.expand_func(derivative).subs(s.log(1 - Z**2), s.log(1 - Z) + s.log(1 + Z))
    assert s.simplify(derivative) == 0
    assert s.limit(result, Z, 0, dir="+") == 0


def test_second_primitive_mobius_equation_exactly():
    hm = HPL((-1,), Z)
    hm1 = HPL((-1, 1), Z)
    expected = (
        hpl_alternating_at_one((-1, 0, 1))
        - s.zeta(2) * hm
        - s.log(2) * hm**2 / 2
        + hm**3 / 6
        - s.log(2) * hm1
        + hm1 * hm
        - HPL((-1, -1, 0), Z)
        - 2 * HPL((-1, -1, 1), Z)
        - HPL((-1, 1, 0), Z)
    )
    assert s.expand(lyndon_reduce(mobius_hpl((-1, 0, 1), Z, **DOMAIN)) - expected) == 0
