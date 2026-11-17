from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl.notation import (
    HPLPM,
    compress_indices,
    expand_indices,
    hpl_pm,
    letters_to_pm,
    pm_to_letters,
)
from hypra.symbolic import HPL


def test_all_short_word_roundtrips():
    for n in range(5):
        for word in product((-1, 0, 1), repeat=n):
            assert expand_indices(compress_indices(word)) == word
            result = {}
            for pm, coefficient in letters_to_pm(word).items():
                for letters, c in pm_to_letters(pm).items():
                    result[letters] = result.get(letters, 0) + coefficient * c
            assert {w: c for w, c in result.items() if c} == {word: 1}


def test_kernel_conventions():
    z = s.Symbol("z")
    assert hpl_pm(("plus",), z) == HPL((1,), z) + HPL((-1,), z)
    assert hpl_pm(("minus",), z) == HPL((1,), z) - HPL((-1,), z)
    assert s.simplify(s.diff(hpl_pm(("plus",), z), z) - 2 / (1 - z**2)) == 0
    assert s.simplify(s.diff(hpl_pm(("minus",), z), z) - 2 * z / (1 - z**2)) == 0
    assert expand_indices((2, -3, 0)) == (0, 1, 0, 0, -1, 0)


@pytest.mark.parametrize("word", [(True,), (0.0,), ("bad",)])
def test_bad_pm_letters(word):
    with pytest.raises((ValueError, TypeError)):
        pm_to_letters(word)


def test_pm_expression_algebra_and_integration():
    from hypra.hpl.algebra import hpl_series, integrate_from_zero, primitive_hpl, shuffle_expand

    z = s.Symbol("z", positive=True)
    plus = HPLPM(("plus",), z)
    assert s.expand(shuffle_expand(plus**2, z) - 2 * hpl_pm(("plus", "plus"), z)) == 0
    integrand = plus / (1 - z**2)
    primitive = primitive_hpl(integrand, z)
    assert s.simplify(s.diff(primitive, z) - integrand.rewrite(HPL)) == 0
    definite = integrate_from_zero(integrand, z)
    assert s.expand(definite - hpl_pm(("plus", "plus"), z) / 2) == 0
    assert hpl_series(plus, z, 5) == 2 * z + 2 * z**3 / 3


def test_pm_expression_log_and_basis_reductions():
    from hypra.hpl.algebra import shuffle_expand
    from hypra.hpl.reduction import extract_hpl_logs, lyndon_reduce

    z = s.Symbol("z", positive=True)
    original = HPLPM(("plus", 0, "minus", 0), z)
    extracted = extract_hpl_logs(original)
    assert not extracted.has(HPLPM)
    assert all(atom.word[0] != 1 and atom.word[-1] != 0 for atom in extracted.atoms(HPL))
    assert shuffle_expand(extracted - original.rewrite(HPL), z) == 0
    reduced = lyndon_reduce(original)
    assert not reduced.has(HPLPM)
    assert shuffle_expand(reduced - original.rewrite(HPL), z) == 0


def test_pm_expression_argument_unification_and_unknown_argument():
    from hypra.errors import UnsupportedExpressionError
    from hypra.hpl.arguments import transform_hpl_arguments
    from hypra.symbolic import reduce_known_hpl

    z = s.Symbol("z", positive=True)
    original = HPLPM(("plus",), -z) + HPLPM((0,), -z)
    unified = transform_hpl_arguments(original, z, domain="0<z<1", boundary="lower")
    expected = -hpl_pm(("plus",), z) + s.log(z) - s.I * s.pi
    assert not unified.has(HPLPM)
    assert s.simplify(reduce_known_hpl(unified - expected)) == 0
    with pytest.raises(UnsupportedExpressionError, match="Unsupported HPL argument"):
        transform_hpl_arguments(HPLPM(("plus",), s.sin(z)), z, domain="0<z<1")


def test_pm_sympy_series_uses_kernel_germ_instead_of_differentiating_word():
    z = s.Symbol("z", positive=True)
    assert HPLPM(("plus",), z).series(z, 0, 5) == 2 * z + 2 * z**3 / 3 + s.Order(z**5, z)
    assert HPLPM((0, "plus"), z).series(z, 0, 5) == (2 * z + 2 * z**3 / 9 + s.Order(z**5, z))
    shifted = HPLPM(("plus", 0), s.Rational(1, 3) + z)
    assert shifted.series(z, 0, 3) == shifted.rewrite(HPL).series(z, 0, 3)


@pytest.mark.parametrize("word", [(0, -1, 1, -1), (0, 1, 1, -1, 1), (0, 1, 1, 0, -1, -1)])
def test_at_i_primitive_evalf_avoids_symbolic_reduction_recursion(word, monkeypatch):
    from hypra.hpl import special
    from hypra.hpl.radial_continuation import continue_hpl

    def no_recursive_reduction(*args, **kwargs):
        pytest.fail("Numerical primitives must bypass the canonical at-i reducer")

    monkeypatch.setattr(special, "hpl_at_i", no_recursive_reduction)
    monkeypatch.setattr(special, "hpl_pm_at_i", no_recursive_reduction)
    expression = HPLPM.from_indices(word, s.I)
    actual = expression.evalf(32)
    assert not actual.has(HPL, HPLPM)

    # Independent continuation integrates the ordinary-letter suffix ODE;
    # it does not use the exact MPL representation behind the primitive.
    ctx = mp.mp.clone()
    ctx.dps = 45
    expected = ctx.mpc(0)
    for ordinary, coefficient in pm_to_letters(expression.word).items():
        value = continue_hpl(
            ordinary,
            s.I,
            domain="principal-cut-plane",
            path="radial",
            dps=32,
        ).value
        expected += int(coefficient) * ctx.mpc(value)
    real, imag = actual.as_real_imag()
    assert abs(ctx.mpc(str(real), str(imag)) - expected) < ctx.mpf("1e-30")


def test_pm_primitive_is_preserved_as_native_expansion_coefficient():
    from hypra import expand_epsilon

    epsilon = s.Symbol("eps")
    primitive = HPLPM((0, "minus", "plus", "minus"), s.I)
    result = expand_epsilon((1 + epsilon) * primitive, epsilon, through=1)
    assert result.coefficient(0) == primitive
    assert result.coefficient(1) == primitive
