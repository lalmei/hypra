"""Independent numeric checks of shuffle and shifted-basepoint reductions."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl.known_functions import reduce_hpl_to_known_functions
from hypra.hpl.notation import HPLPM, expand_indices
from hypra.hpl.numeric import evaluate_hpl
from hypra.nielsen import NielsenPolylog, nielsen_to_hpl
from hypra.symbolic import HPL

Z = s.Symbol("z", positive=True)


def numerical(expression, point, digits=32):
    expression = nielsen_to_hpl(expression).subs(Z, point)
    # For Nielsen words there are no trailing zeros: negation therefore has
    # no log(z) branch term and follows directly by changing t -> -t.
    expression = expression.xreplace(
        {
            atom: (-1) ** sum(letter == 1 for letter in atom.word)
            * HPL(tuple(-int(letter) for letter in atom.word), -atom.argument)
            for atom in expression.atoms(HPL)
            if atom.argument.is_negative
        }
    )
    return expression.evalf(digits)


@pytest.mark.parametrize(
    "word", [word for word in product((-1, 0, 1), repeat=4) if len(set(word)) <= 2], ids=str
)
def test_every_two_letter_quartic_against_interior_series(word):
    result = reduce_hpl_to_known_functions(HPL(word, Z), Z, domain="0<z<1")
    assert not result.has(HPL)
    context = mp.mp.clone()
    context.dps = 35
    for point in (s.Rational(1, 3), s.Rational(3, 4)):
        actual = context.mpf(str(numerical(result, point)))
        reference = evaluate_hpl(word, point, dps=28).value
        assert abs(actual - reference) < context.mpf("1e-27")


@pytest.mark.parametrize(
    "indices",
    [
        (1, 3),
        (1, 1, 2),
        (1, 2, 1),
        (2, 2),
        (3, 1),
        (2, 1, 1),
        (1, 1, 1, 1),
        (-1, 1, 1, 1),
        (-1, 1, -1, -1),
        (-3, -1),
        (-2, -1, -1),
        (4, 1),
        (3, 1, 1),
        (5, 1, 1, 1),
    ],
)
def test_pinned_upstream_higher_weight_known_function_patterns(indices):
    value = HPL(expand_indices(indices), Z)
    assert not reduce_hpl_to_known_functions(value, Z, domain="0<z<1").has(HPL)


def test_signed_nielsen_family_has_the_correct_sign():
    assert reduce_hpl_to_known_functions(HPL((0, 0, -1, -1, -1), Z), Z, domain="0<z<1") == (
        -NielsenPolylog(2, 3, -Z)
    )


def test_plus_minus_known_function_dispatch():
    assert (
        s.simplify(
            reduce_hpl_to_known_functions(HPLPM(("plus",), Z), Z, domain="0<z<1")
            - s.log(1 + Z)
            + s.log(1 - Z)
        )
        == 0
    )


def test_unknown_mixed_words_are_retained_and_domains_are_checked():
    value = HPL((-1, 0, 1, 0), Z)
    assert reduce_hpl_to_known_functions(value, Z, domain="0<z<1") == value
    with pytest.raises(UnsupportedExpressionError, match="domain"):
        reduce_hpl_to_known_functions(value, Z, domain="complex")
    with pytest.raises(UnsupportedExpressionError, match="common argument"):
        reduce_hpl_to_known_functions(HPL((-1, 0, 1, 0), 1 - Z), Z, domain="0<z<1")
