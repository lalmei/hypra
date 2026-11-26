import pytest
import sympy as s

from hypra import HPL, UnsupportedExpressionError, hpl_at_minus_one


def test_finite_minus_one_endpoints_and_sympy_evalf():
    for word, expected in [
        ((), 1),
        ((1,), -s.log(2)),
        ((0, -1), -s.zeta(2)),
        ((0, 1, -1), 13 * s.zeta(3) / 8 - s.pi**2 * s.log(2) / 4),
        ((0, 0), -(s.pi**2) / 2),
    ]:
        assert hpl_at_minus_one(word) == expected
        assert hpl_at_minus_one(word, boundary="upper") == hpl_at_minus_one(word, boundary="lower")
        actual = HPL(word, -1).evalf(40)
        assert abs(actual - s.N(expected, 50)) < s.Rational(1, 10**39)


@pytest.mark.parametrize(
    "word,upper",
    [
        ((0,), s.I * s.pi),
        ((0,) * 23, (s.I * s.pi) ** 23 / s.factorial(23)),
        ((1, 0), s.zeta(2) / 2 - s.I * s.pi * s.log(2)),
        ((0, 1, 0), 3 * s.zeta(3) / 2 - s.I * s.pi * s.zeta(2) / 2),
    ],
)
def test_trailing_zeros_require_a_side(word, upper):
    assert s.expand(hpl_at_minus_one(word, boundary="upper") - upper) == 0
    assert s.expand(hpl_at_minus_one(word, boundary="lower") - s.conjugate(upper)) == 0
    with pytest.raises(UnsupportedExpressionError, match="boundary"):
        hpl_at_minus_one(word)
    assert HPL(word, -1).evalf(30) == HPL(word, -1)


@pytest.mark.parametrize("word", [(-1,), (-1, 0), (-1, 1, 0), (-1, -1, 0, 0)])
def test_divergent_minus_one_endpoints(word):
    for side in ("upper", "lower"):
        with pytest.raises(UnsupportedExpressionError, match="Divergent|No finite"):
            hpl_at_minus_one(word, boundary=side)
    assert HPL(word, -1).evalf() == HPL(word, -1)


def test_minus_one_validates_input_and_weight():
    with pytest.raises(ValueError, match="boundary"):
        hpl_at_minus_one((1,), boundary="both")
    with pytest.raises(UnsupportedExpressionError):
        hpl_at_minus_one((1.0,))
    with pytest.raises(UnsupportedExpressionError, match="weight"):
        hpl_at_minus_one((0,) * 25)
