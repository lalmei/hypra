"""Canonical at-i identities checked without imported reduction coefficients."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.hpl.at_i_reduction import _PRIMITIVES, _mobius_terms, at_i_weight_reductions
from hypra.hpl.notation import HPLPM, _pm_indices, pm_to_letters
from hypra.hpl.radial_continuation import continue_hpl
from hypra.hpl.special import _ordinary_at_i_raw, _reduce_at_i_constants, hpl_at_i, hpl_pm_at_i
from hypra.symbolic import HPL


@pytest.mark.parametrize("weight", range(7))
def test_complete_weight_space_uses_only_declared_canonical_primitives(weight):
    table = at_i_weight_reductions(weight)
    assert len(table) == 3**weight
    expected = {
        HPLPM.from_indices(word, s.I)
        for degree, words in _PRIMITIVES.items()
        if degree <= weight
        for word in words
    }
    actual = set().union(*(value.atoms(HPLPM) for value in table.values()))
    assert actual == expected
    assert all(not value.has(HPL, s.Float) for value in table.values())
    for word in _PRIMITIVES.get(weight, ()):
        assert table[word] == HPLPM.from_indices(word, s.I)


def test_reduction_table_is_immutable():
    table = at_i_weight_reductions(4)
    word = (1, 0, 0, 1)
    expected = table[word]
    with pytest.raises(TypeError):
        table[word] = s.S.Zero
    assert table[word] == expected


def test_quartic_shuffle_identity_closes_to_named_constants():
    # Integrating H_+ H_00+ by parts, then using the shuffle square of H_0+,
    # gives H_+00+=H_+H_00+-(H_0+)²/2 for any nonsingular argument.
    value = hpl_pm_at_i(("plus", 0, 0, "plus"))
    assert value == 2 * s.Catalan**2 - s.pi**4 / 32
    assert (
        s.expand(
            value
            - hpl_pm_at_i(("plus",)) * hpl_pm_at_i((0, 0, "plus"))
            + hpl_pm_at_i((0, "plus")) ** 2 / 2
        )
        == 0
    )


def test_cubic_conjugate_normalization_matches_independently_derived_formulas():
    pair = s.polylog(3, (1 + s.I) / 2) + s.polylog(3, (1 - s.I) / 2)
    assert _reduce_at_i_constants(pair) == (
        -5 * s.pi**2 * s.log(2) / 96 + s.log(2) ** 3 / 24 + 35 * s.zeta(3) / 32
    )
    for word in product((-1, 0, 1), repeat=3):
        assert s.expand(hpl_at_i(word) - _ordinary_at_i_raw(word)) == 0


@pytest.mark.parametrize("weight", [4, 5, 6])
def test_repeated_letters_and_classical_pm_values(weight):
    assert hpl_pm_at_i(("plus",) * weight) == (s.I * s.pi / 2) ** weight / s.factorial(weight)
    assert hpl_pm_at_i(("minus",) * weight) == (-s.log(2)) ** weight / s.factorial(weight)
    assert (
        s.expand(
            hpl_pm_at_i((0,) * (weight - 1) + ("plus",))
            - hpl_at_i((0,) * (weight - 1) + (1,))
            - hpl_at_i((0,) * (weight - 1) + (-1,))
        )
        == 0
    )


@pytest.mark.parametrize(
    "word",
    [
        (-1, 0),
        (1, -1, 0),
        (0, 1, -1, -1),
        (-1, 0, 1, 0),
        (0, -1, 1, -1, -1),
    ],
)
def test_mobius_integration_constants_against_independent_interior_values(word):
    context = mp.mp.clone()
    context.dps = 36
    for point in (s.Rational(1, 3), s.Rational(3, 4)):
        transformed = s.Add(
            *(
                coefficient * (HPLPM.from_indices(target, point) if target else s.S.One)
                for target, coefficient in _mobius_terms(word)
            )
        )
        original = HPLPM.from_indices(word, (1 - point) / (1 + point))
        gap = context.mpf(str((transformed - original).evalf(32)))
        assert abs(gap) < context.mpf("1e-28")


@pytest.mark.parametrize(
    "word",
    [
        (0, -1, 1, -1),
        (1, -1, 0, 1),
        (0, -1, 1, -1, -1),
        (1, 0, 0, 1, 0),
        (0, 1, -1, -1, 0, -1),
        (0, 1, 1, 0, 0, 1),
    ],
)
def test_retained_and_reduced_values_against_independent_radial_continuation(word):
    context = mp.mp.clone()
    context.dps = 37
    pm_word = HPLPM.from_indices(word, s.I).word
    source = hpl_pm_at_i(pm_word)
    real, imag = source.evalf(33).as_real_imag()
    actual = context.mpc(str(real), str(imag))
    reference = context.zero
    for letters, coefficient in pm_to_letters(pm_word).items():
        reference += (
            int(coefficient)
            * continue_hpl(letters, s.I, domain="principal-cut-plane", path="radial", dps=28).value
        )
    assert abs(actual - reference) < context.mpf("1e-26")


def test_ordinary_and_pm_quartic_functional_identity():
    # The linear kernel relation holds before integration and at every point.
    # Checking its endpoint after reduction also exercises cancellations of
    # the retained primitive, rather than just numerical agreement.
    total = s.Add(*((-1) ** word.count(-1) * hpl_at_i(word) for word in product((-1, 1), repeat=4)))
    assert s.expand(total - hpl_pm_at_i(("minus",) * 4)) == 0
    assert _pm_indices(HPLPM((0, "minus", "plus", "minus"), s.I).word) == (0, -1, 1, -1)
