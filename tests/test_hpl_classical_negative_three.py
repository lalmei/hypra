"""Real classical reduction of all weight-three {0,-1} words."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl_classical import classical_hpl
from hypra.hpl_numeric import evaluate_hpl
from hypra.symbolic import HPL

Z = s.Symbol("z", positive=True)
DOMAIN = {"domain": "0<z<1"}


@pytest.mark.parametrize("word", list(product((0, -1), repeat=3)))
def test_negative_cubic_classical_derivative_and_values(word):
    expression = classical_hpl(HPL(word, Z), Z, **DOMAIN)
    assert not expression.has(HPL, s.I)
    kernel = 1 / Z if word[0] == 0 else 1 / (1 + Z)
    tail = classical_hpl(HPL(word[1:], Z), Z, **DOMAIN)
    residual = s.expand_func(s.diff(expression, Z) - kernel * tail)
    # All logarithm arguments here are provably positive for positive Z.
    assert s.simplify(s.expand_log(residual, force=False)) == 0
    context = mp.mp.clone()
    context.dps = 40
    for point in (s.Rational(1, 3), s.Rational(3, 4)):
        actual = context.mpf(str(s.N(expression.subs(Z, point), 38)))
        reference = evaluate_hpl(word, point, dps=30).value
        assert abs(actual - reference) < context.mpf("1e-29")


def test_negative_new_formula_origin_and_real_arguments():
    expression = classical_hpl(HPL((0, -1, -1), Z), Z, **DOMAIN)
    assert s.limit(expression, Z, 0, dir="+") == 0
    assert {p.args[1] for p in expression.atoms(s.polylog)} == {1 / (1 + Z)}


def test_true_ternary_cubic_stays_lossless():
    atom = HPL((1, 0, -1, 0), Z)
    assert classical_hpl(atom, Z, **DOMAIN) == atom
