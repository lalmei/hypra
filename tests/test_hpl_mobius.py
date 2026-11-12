"""Pullback differential equations and independent elementary special cases."""

from itertools import product

import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import shuffle_expand
from hypra.hpl_mobius import mobius_hpl
from hypra.symbolic import HPL, reduce_known_hpl

Z = s.Symbol("z", positive=True)
DOMAIN = {"domain": "0<z<1"}


@pytest.mark.parametrize("word", [w for n in range(1, 4) for w in product((-1, 0, 1), repeat=n)])
def test_mobius_differential_pullback(word):
    kernel = {0: -1 / (1 - Z) - 1 / (1 + Z), 1: -1 / Z + 1 / (1 + Z), -1: -1 / (1 + Z)}[word[0]]
    transformed = mobius_hpl(word, Z, **DOMAIN)
    residual = s.diff(transformed, Z) - kernel * mobius_hpl(word[1:], Z, **DOMAIN)
    assert s.simplify(shuffle_expand(residual, Z)) == 0
    assert all(atom.argument == Z for atom in transformed.atoms(HPL))


@pytest.mark.parametrize(
    "letter,expected",
    [
        (0, s.log(1 - Z) - s.log(1 + Z)),
        (1, -s.log(2) - s.log(Z) + s.log(1 + Z)),
        (-1, s.log(2) - s.log(1 + Z)),
    ],
)
def test_mobius_weight_one_constants(letter, expected):
    assert s.simplify(reduce_known_hpl(mobius_hpl((letter,), Z, **DOMAIN)) - expected) == 0


def test_repeated_letters_and_domain_limits():
    for letter in (-1, 0, 1):
        one = mobius_hpl((letter,), Z, **DOMAIN)
        residual = mobius_hpl((letter,) * 3, Z, **DOMAIN) - one**3 / 6
        assert s.simplify(shuffle_expand(residual, Z)) == 0
    with pytest.raises(UnsupportedExpressionError):
        mobius_hpl((0,), Z)
    with pytest.raises(UnsupportedExpressionError):
        mobius_hpl((0,) * 7, Z, **DOMAIN)


@pytest.mark.parametrize("word", [(0, -1, 1, 0), (1, -1, 0, 1, 0), (-1, 0, -1, 1, 1)])
def test_higher_weight_mobius_pullback(word):
    test_mobius_differential_pullback(word)
