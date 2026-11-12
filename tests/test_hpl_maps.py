import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import hpl_series, shuffle_expand
from hypra.hpl_maps import fractional_hpl, square_hpl
from hypra.symbolic import HPL, reduce_known_hpl

Z = s.Symbol("z", positive=True)
DOMAIN = {"domain": "0<z<1"}
WORDS = [(), (0,), (1,), (0, 0), (0, 1), (1, 0), (1, 1), (1, 0, 1), (0, 1, 0), (0, 0, 1, 1)]


@pytest.mark.parametrize("word", WORDS)
def test_square_map_exact_local_series(word):
    actual = hpl_series(square_hpl(word, Z, **DOMAIN), Z, 8)
    original = hpl_series(HPL(word, Z), Z, 4).subs(Z, Z**2)
    original = s.expand_log(original, force=False)
    assert s.simplify(actual - original) == 0


@pytest.mark.parametrize("word", WORDS[1:])
@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
def test_fractional_map_chain_rule_and_origin(word, boundary, sign):
    result = fractional_hpl(word, Z, boundary=boundary, **DOMAIN)
    tail = fractional_hpl(word[1:], Z, boundary=boundary, **DOMAIN)
    kernel = 1 / Z + 1 / (1 - Z) if word[0] == 0 else -1 / (1 - Z)
    assert s.simplify(shuffle_expand(s.diff(result, Z) - kernel * tail, Z)) == 0
    if word[-1] != 0:
        assert s.limit(hpl_series(result, Z, 2), Z, 0) == 0
    elif all(letter == 0 for letter in word):
        expected = (s.log(Z) - s.log(1 - Z) + sign * s.I * s.pi) ** len(word) / s.factorial(
            len(word)
        )
        assert s.simplify(reduce_known_hpl(result) - expected) == 0


def test_fractional_dilog_identity_and_monodromy():
    expected = -HPL((0, 1), Z) - HPL((1,), Z) ** 2 / 2
    for boundary in ("upper", "lower"):
        assert (
            shuffle_expand(fractional_hpl((0, 1), Z, boundary=boundary, **DOMAIN) - expected, Z)
            == 0
        )
    jump = fractional_hpl((1, 0), Z, boundary="upper", **DOMAIN) - fractional_hpl(
        (1, 0), Z, boundary="lower", **DOMAIN
    )
    assert s.expand(jump + 2 * s.I * s.pi * HPL((1,), Z)) == 0


def test_map_rejections():
    for func in (square_hpl, lambda w, z, **kw: fractional_hpl(w, z, boundary="upper", **kw)):
        with pytest.raises(UnsupportedExpressionError):
            func((-1,), Z, **DOMAIN)
        with pytest.raises(UnsupportedExpressionError):
            func((0,) * 11, Z, **DOMAIN)
        with pytest.raises(UnsupportedExpressionError):
            func((0,), Z)
    with pytest.raises(UnsupportedExpressionError):
        fractional_hpl((1,), Z, boundary="unspecified", **DOMAIN)
