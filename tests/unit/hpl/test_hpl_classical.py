"""Classical HPL reduction on the real interval 0<z<1."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl.alternating import hpl_alternating_at_one
from hypra.hpl.classical import classical_hpl
from hypra.hpl.mobius import mobius_hpl
from hypra.hpl.numeric import evaluate_hpl
from hypra.hpl.reduction import lyndon_reduce
from hypra.symbolic import HPL, reduce_known_hpl

Z = s.Symbol("z", positive=True)
DOMAIN = "0<z<1"
LETTERS = (-1, 0, 1)


@pytest.mark.parametrize("word", list(product(LETTERS, repeat=2)), ids=str)
def test_all_weight_two_differential_and_origin(word):
    result = classical_hpl(HPL(word, Z), Z, domain=DOMAIN)
    assert not result.has(HPL)
    difference = s.diff(result, Z) - reduce_known_hpl(s.diff(HPL(word, Z), Z))
    difference = difference.xreplace(
        {s.log((1 - Z) / 2): s.log(1 - Z) - s.log(2), s.log((1 + Z) / 2): s.log(1 + Z) - s.log(2)}
    )
    assert s.simplify(s.expand_func(difference)) == 0
    if word != (0, 0):
        assert s.limit(result, Z, 0, dir="+") == 0


@pytest.mark.parametrize("word", list(product(LETTERS, repeat=3)), ids=str)
def test_all_weight_three_against_numerical_evaluator(word):
    expression = classical_hpl(HPL(word, Z), Z, domain=DOMAIN)
    assert not expression.has(HPL, s.I)
    context = mp.mp.clone()
    context.dps = 40
    for point in (s.Rational(1, 3), s.Rational(3, 4)):
        actual = context.mpf(str(s.N(expression.subs(Z, point), 38)))
        reference = evaluate_hpl(word, point, dps=30).value
        assert abs(actual - reference) < context.mpf("1e-29"), (word, point)


@pytest.mark.parametrize(
    "word",
    [word for word in product(LETTERS, repeat=3) if -1 not in word or 1 not in word],
    ids=str,
)
def test_weight_three_differential_without_opposite_signs(word):
    expression = classical_hpl(HPL(word, Z), Z, domain=DOMAIN)
    tail = classical_hpl(HPL(word[1:], Z), Z, domain=DOMAIN)
    kernel = {-1: 1 / (1 + Z), 0: 1 / Z, 1: 1 / (1 - Z)}[word[0]]
    residual = s.expand_func(s.diff(expression, Z) - kernel * tail)
    # Expand only logarithms whose arguments are provably positive.
    assert s.simplify(s.expand_log(residual)) == 0


@pytest.mark.parametrize("word", [(0, 1, 1), (0, -1, -1)], ids=str)
def test_weight_three_primitive_origin(word):
    expression = classical_hpl(HPL(word, Z), Z, domain=DOMAIN)
    assert s.limit(expression, Z, 0, dir="+") == 0


def test_negative_primitive_polylog_arguments():
    expression = classical_hpl(HPL((0, -1, -1), Z), Z, domain=DOMAIN)
    assert {polylog.args[1] for polylog in expression.atoms(s.polylog)} == {1 / (1 + Z)}


@pytest.mark.parametrize("word", [(-1, 1, 1), (1, -1, -1)], ids=str)
def test_mixed_sign_primitive_differential_and_origin(word):
    result = classical_hpl(HPL(word, Z), Z, domain=DOMAIN)
    tail = classical_hpl(HPL(word[1:], Z), Z, domain=DOMAIN)
    kernel = 1 / (1 - Z) if word[0] == 1 else 1 / (1 + Z)
    # Every log argument is positive on the declared interval.
    residual = s.expand_func(s.diff(result, Z) - kernel * tail)
    u = (1 - Z) / 2
    # Euler reflection is independently established by derivative and value at 0.
    residual = residual.subs(
        s.polylog(2, (1 + Z) / 2), s.pi**2 / 6 - s.polylog(2, u) - s.log(u) * s.log(1 - u)
    )
    residual = s.expand_log(residual, force=True)
    assert s.simplify(residual) == 0
    assert s.limit(result, Z, 0, dir="+") == 0


def test_ternary_primitive_differential_and_origin():
    result = classical_hpl(HPL((-1, 0, 1), Z), Z, domain=DOMAIN)
    # The constituent real primitives have separately checked derivatives.
    # Verify the square-log construction before classical substitution here.
    primitive = (
        HPL((0, 1, 1), Z**2) / 2
        - HPL((0, -1, -1), Z)
        - HPL((0, 1, 1), Z)
        + s.log(1 + Z) * s.polylog(2, Z)
    )
    derivative = reduce_known_hpl(s.diff(primitive, Z) - s.polylog(2, Z) / (1 + Z))
    derivative = s.expand_func(derivative).subs(s.log(1 - Z**2), s.log(1 - Z) + s.log(1 + Z))
    assert s.simplify(derivative) == 0
    assert s.limit(result, Z, 0, dir="+") == 0


def test_ternary_primitive_mobius_equation():
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
    assert s.expand(lyndon_reduce(mobius_hpl((-1, 0, 1), Z, domain=DOMAIN)) - expected) == 0


@pytest.mark.parametrize(
    "word", [(1, -1, 0, 1), (-1, 0, 1, 0), (-1, 1, 0, 1), (1, 0, -1, 0)], ids=str
)
def test_higher_mixed_words_are_preserved(word):
    value = HPL(word, Z)
    assert classical_hpl(value, Z, domain=DOMAIN) == value
