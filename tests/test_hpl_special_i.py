from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_special import hpl_at_i
from hypra.symbolic import HPL


@pytest.mark.parametrize("word", list(product((-1, 0, 1), repeat=2)))
def test_all_weight_two_at_i_against_independent_path_integral(word):
    context = mp.mp.clone()
    context.dps = 40

    def integrand(t):
        z = context.j * t
        tail = {0: context.log(z), 1: -context.log(1 - z), -1: context.log(1 + z)}[word[1]]
        kernel = {0: 1 / z, 1: 1 / (1 - z), -1: 1 / (1 + z)}[word[0]]
        return context.j * kernel * tail

    if word == (0, 0):
        reference = (
            -(context.pi**2) / 8
        )  # logarithmic tangential origin, not a convergent path integral
    else:
        reference = context.quad(integrand, [0, 1])
    actual = hpl_at_i(word)
    assert not actual.has(HPL, s.Float)
    real, imag = actual.as_real_imag()
    numerical = context.mpc(str(s.N(real, 38)), str(s.N(imag, 38)))
    assert abs(numerical - reference) < context.mpf("1e-35")


def test_principal_logs_and_catalan_constants():
    assert hpl_at_i(()) == 1
    assert hpl_at_i((0,)) == s.I * s.pi / 2
    assert s.simplify(hpl_at_i((1,)) + s.log(2) / 2 - s.I * s.pi / 4) == 0
    assert s.simplify(hpl_at_i((-1,)) - s.log(2) / 2 - s.I * s.pi / 4) == 0
    assert hpl_at_i((0, 1)) == -(s.pi**2) / 48 + s.I * s.Catalan
    assert hpl_at_i((0, -1)) == s.pi**2 / 48 + s.I * s.Catalan


def test_simple_higher_words_and_rejections():
    assert hpl_at_i((0, 0, 1)) == s.polylog(3, s.I)
    assert hpl_at_i((0, 0, -1)) == -s.polylog(3, -s.I)
    assert s.simplify(hpl_at_i((1,) * 4) - hpl_at_i((1,)) ** 4 / 24) == 0
    with pytest.raises(UnsupportedExpressionError):
        hpl_at_i((1, 0, -1))
    with pytest.raises(UnsupportedExpressionError):
        hpl_at_i((0,) * 25)
