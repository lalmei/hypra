import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_continuation import analytic_continue_hpl
from hypra.symbolic import HPL, reduce_known_hpl

X = s.Symbol("x", real=True)
REGIONS = {
    "0<x<1": s.Rational(1, 2),
    "-1<x<0": -s.Rational(1, 2),
    "x>1": s.Integer(2),
    "x<-1": -s.Integer(2),
}


@pytest.mark.parametrize("region", REGIONS)
@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
@pytest.mark.parametrize("letter", [-1, 0, 1])
def test_weight_one_boundary_signs(region, boundary, sign, letter):
    value = REGIONS[region]
    result = reduce_known_hpl(
        analytic_continue_hpl((letter,), X, region=region, boundary=boundary)
    ).subs(X, value)
    argument = {0: value, 1: 1 - value, -1: 1 + value}[letter]
    coefficient = -1 if letter == 1 else 1
    imaginary = 0
    if argument < 0:
        imaginary = sign * s.I * s.pi * (-1 if letter == 1 else 1)
    expected = coefficient * (s.log(abs(argument)) + imaginary)
    assert s.simplify(result - expected) == 0


@pytest.mark.parametrize("boundary", ["upper", "lower"])
def test_negative_far_dilog_is_real_and_correct(boundary):
    result = reduce_known_hpl(
        analytic_continue_hpl((0, 1), X, region="x<-1", boundary=boundary)
    ).subs(X, -2)
    # Independent classical inversion Li2(-2)+Li2(-1/2)=-pi²/6-log²2/2.
    expected = -s.polylog(2, -s.Rational(1, 2)) - s.pi**2 / 6 - s.log(2) ** 2 / 2
    # Remaining HPLs are at positive 1/2 and reduce by shuffle extraction.
    from hypra.hpl_reduce import extract_hpl_logs

    result = extract_hpl_logs(result)
    assert abs(s.N(result - expected, 35)) < s.Float("1e-30")


def test_region_contract_and_rejections():
    assert analytic_continue_hpl((0, 1), X, region="0<x<1", boundary="upper") == HPL((0, 1), X)
    for region, boundary in [("unknown", "upper"), ("x>1", "unknown")]:
        with pytest.raises(UnsupportedExpressionError):
            analytic_continue_hpl((0,), X, region=region, boundary=boundary)
    with pytest.raises(UnsupportedExpressionError):
        analytic_continue_hpl((0,) * 7, X, region="x>1", boundary="upper")
