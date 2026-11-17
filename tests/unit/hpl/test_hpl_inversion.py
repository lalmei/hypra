from itertools import product

import pytest
import sympy as s

from hypra.hpl.algebra import shuffle_expand
from hypra.hpl.inversion import invert_hpl
from hypra.symbolic import HPL, reduce_known_hpl

Z = s.Symbol("z", positive=True)
DOMAIN = {"domain": "0<z<1"}


@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
@pytest.mark.parametrize("letter", [-1, 0, 1])
def test_weight_one_inversion_boundaries(boundary, sign, letter):
    expected = {
        -1: s.log(1 + Z) - s.log(Z),
        0: -s.log(Z),
        1: s.log(Z) - s.log(1 - Z) + sign * s.I * s.pi,
    }[letter]
    actual = reduce_known_hpl(invert_hpl((letter,), Z, boundary=boundary, **DOMAIN))
    assert s.simplify(actual - expected) == 0


@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
def test_dilogarithm_inversion(boundary, sign):
    expected = -s.polylog(2, Z) + s.pi**2 / 3 - s.log(Z) ** 2 / 2 - sign * s.I * s.pi * s.log(Z)
    actual = invert_hpl((0, 1), Z, boundary=boundary, **DOMAIN)
    assert s.simplify(shuffle_expand(actual - expected, Z)) == 0


@pytest.mark.parametrize("word", list(product((-1, 0, 1), repeat=3)))
@pytest.mark.parametrize("boundary", ["upper", "lower"])
def test_inversion_chain_rule(word, boundary):
    kernel = {0: -1 / Z, 1: 1 / Z + 1 / (1 - Z), -1: -1 / Z + 1 / (1 + Z)}[word[0]]
    actual = invert_hpl(word, Z, boundary=boundary, **DOMAIN)
    tail = invert_hpl(word[1:], Z, boundary=boundary, **DOMAIN)
    assert s.simplify(shuffle_expand(s.diff(actual, Z) - kernel * tail, Z)) == 0


def test_inversion_input_rejections():
    from hypra.errors import UnsupportedExpressionError

    with pytest.raises(UnsupportedExpressionError):
        invert_hpl((0,), Z, boundary="upper")
    with pytest.raises(UnsupportedExpressionError):
        invert_hpl((0,), Z, boundary="unknown", **DOMAIN)
    with pytest.raises(UnsupportedExpressionError):
        invert_hpl((0,) * 8 + (1,), Z, boundary="lower", **DOMAIN)


@pytest.mark.parametrize("letter", [-1, 0, 1])
@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
def test_repeated_letters_through_weight_twenty_four(letter, boundary, sign):
    primitive = {
        -1: s.log(1 + Z) - s.log(Z),
        0: -s.log(Z),
        1: s.log(Z) - s.log(1 - Z) + sign * s.I * s.pi,
    }[letter]
    actual = reduce_known_hpl(invert_hpl((letter,) * 24, Z, boundary=boundary, **DOMAIN))
    assert actual == primitive**24 / s.factorial(24)


@pytest.mark.parametrize("degree", range(2, 9))
@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
def test_classical_inversion_through_weight_eight(degree, boundary, sign):
    # Independent classical inversion in terms of Bernoulli polynomials.
    word = (0,) * (degree - 1) + (1,)
    phase = s.Rational(1 - sign, 2) - s.log(Z) / (2 * s.pi * s.I)
    expected = -((-1) ** degree) * s.polylog(degree, Z) - (2 * s.pi * s.I) ** degree / s.factorial(
        degree
    ) * s.bernoulli(degree, phase)
    actual = reduce_known_hpl(invert_hpl(word, Z, boundary=boundary, **DOMAIN))
    assert s.simplify(actual - expected) == 0


@pytest.mark.parametrize("word", list(product((-1, 0, 1), repeat=3)))
@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
def test_inversion_endpoint_matches_branch_log_shift(word, boundary, sign):
    from hypra.hpl.regularization import (
        hpl_endpoint_polynomial,
        hpl_expression_endpoint_polynomial,
    )

    T = s.Symbol("T")
    actual = invert_hpl(word, Z, boundary=boundary, **DOMAIN)
    endpoint = hpl_expression_endpoint_polynomial(actual, Z, T)
    expected = hpl_endpoint_polynomial(word, T).subs(T, T + sign * s.I * s.pi)
    assert s.simplify(endpoint - expected) == 0


def test_weight_eight_mixed_inversion_chain_rule():
    word = (-1, 0, 1, 0, -1, 0, 1, 0)
    actual = invert_hpl(word, Z, boundary="upper", **DOMAIN)
    tail = invert_hpl(word[1:], Z, boundary="upper", **DOMAIN)
    residual = s.expand(s.diff(actual, Z) - (-1 / Z + 1 / (1 + Z)) * tail)
    # Compare coefficients in the formal HPL basis; cancellation of named
    # endpoint constants must not depend on a numerical sample.
    for coefficient in s.collect(residual, list(residual.atoms(HPL)), evaluate=False).values():
        assert s.simplify(coefficient) == 0
