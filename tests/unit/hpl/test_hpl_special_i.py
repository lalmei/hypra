from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl.special import hpl_at_i
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
    assert hpl_at_i((0, 0, 1)) == -3 * s.zeta(3) / 32 + s.I * s.pi**3 / 32
    assert hpl_at_i((0, 0, -1)) == 3 * s.zeta(3) / 32 + s.I * s.pi**3 / 32
    assert s.simplify(hpl_at_i((1,) * 4) - hpl_at_i((1,)) ** 4 / 24) == 0
    with pytest.raises(UnsupportedExpressionError):
        hpl_at_i((1, 0, -1, 0, 1, 0, 1))
    with pytest.raises(UnsupportedExpressionError):
        hpl_at_i((0,) * 25)


def test_mixed_quartic_uses_the_canonical_at_i_primitive():
    from hypra.hpl.multiple_polylog import MultiplePolylog
    from hypra.hpl.notation import HPLPM

    value = hpl_at_i((1, -1, 0, 1))
    assert value.has(HPLPM)
    assert value.atoms(HPLPM) == {HPLPM((0, "minus", "plus", "minus"), s.I)}
    assert not value.has(HPL, MultiplePolylog, s.Float)
    assert not value.free_symbols


@pytest.mark.parametrize("word", list(product((-1, 0, 1), repeat=3)))
def test_cubic_at_i_named_constants_against_radial_continuation(word):
    from hypra.hpl.multiple_polylog import MultiplePolylog
    from hypra.hpl.radial_continuation import continue_hpl

    result = hpl_at_i(word)
    assert not result.has(HPL, MultiplePolylog, s.Float)
    assert all(
        atom.args[0] == 3 and atom.args[1] in ((1 + s.I) / 2, (1 - s.I) / 2)
        for atom in result.atoms(s.polylog)
    )
    context = mp.mp.clone()
    context.dps = 40
    real, imag = result.evalf(38).as_real_imag()
    actual = context.mpc(str(real), str(imag))
    reference = continue_hpl(word, s.I, domain="principal-cut-plane", path="radial", dps=32).value
    assert abs(actual - reference) < context.mpf("1e-31")


@pytest.mark.parametrize("order", range(2, 13))
def test_classical_quarter_period_constants(order):
    context = mp.mp.clone()
    context.dps = 40
    for sign in (-1, 1):
        value = hpl_at_i((0,) * (order - 1) + (sign,))
        assert not value.has(s.polylog)
        real, imag = value.evalf(38).as_real_imag()
        actual = context.mpc(str(real), str(imag))
        reference = sign * context.polylog(order, sign * context.j)
        assert abs(actual - reference) < context.mpf("1e-35")
