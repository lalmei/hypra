import pytest
import sympy as s

from hypra.hpl_algebra import shuffle_expand
from hypra.hpl_inversion import invert_hpl
from hypra.symbolic import reduce_known_hpl

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


@pytest.mark.parametrize("word", [(-1, 0), (1, -1), (1, 0, 1), (0, -1, 1)])
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
        invert_hpl((0,) * 7, Z, boundary="lower", **DOMAIN)
