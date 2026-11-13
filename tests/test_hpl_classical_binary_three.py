"""Independent derivative/origin and numerical checks of binary weight three."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl_classical import classical_hpl
from hypra.hpl_numeric import evaluate_hpl
from hypra.symbolic import HPL

Z = s.Symbol("z", positive=True)
DOMAIN = {"domain": "0<z<1"}


@pytest.mark.parametrize("word", list(product((0, 1), repeat=3)))
def test_binary_cubic_derivative_and_numeric(word):
    expression = classical_hpl(HPL(word, Z), Z, **DOMAIN)
    assert not expression.has(HPL)
    tail = classical_hpl(HPL(word[1:], Z), Z, **DOMAIN)
    kernel = 1 / Z if word[0] == 0 else 1 / (1 - Z)
    assert s.simplify(s.expand_func(s.diff(expression, Z) - kernel * tail)) == 0
    context = mp.mp.clone()
    context.dps = 40
    for point in (s.Rational(1, 3), s.Rational(3, 4)):
        reference = evaluate_hpl(word, point, dps=30).value
        actual = context.mpf(str(s.N(expression.subs(Z, point), 38)))
        assert abs(actual - reference) < context.mpf("1e-29")


def test_new_irreducible_formula_origin_constant():
    expression = classical_hpl(HPL((0, 1, 1), Z), Z, **DOMAIN)
    assert s.limit(expression, Z, 0, dir="+") == 0


def test_higher_mixed_words_remain_lossless():
    atom = HPL((-1, 0, 1, 0), Z)
    assert classical_hpl(atom, Z, **DOMAIN) == atom
