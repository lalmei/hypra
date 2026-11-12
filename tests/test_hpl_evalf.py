"""Opt-in SymPy numeric evaluation without implicit complex continuation."""

import pytest
import sympy as s

from hypra import HPL


@pytest.mark.parametrize('word, expected', [
    ((0,), s.log(s.Rational(1, 3))),
    ((0, 1), s.polylog(2, s.Rational(1, 3))),
    ((1, 0), -s.log(s.Rational(1, 3)) * s.log(s.Rational(2, 3))
     - s.polylog(2, s.Rational(1, 3))),
    ((0, 0, -1), -s.polylog(3, -s.Rational(1, 3))),
])
def test_explicit_evalf_preserves_requested_precision(word, expected):
    expr = HPL(word, s.Rational(1, 3))
    assert isinstance(expr, HPL)
    result = expr.evalf(50)
    assert isinstance(result, s.Float)
    assert abs(result - expected.evalf(70)) < s.Float('1e-49')


@pytest.mark.parametrize('argument', [s.Symbol('z'), s.I / 3, -s.Rational(1, 3), 0, 1, 2])
def test_evalf_unknown_or_outside_real_domain_stays_symbolic(argument):
    assert HPL((1, 0), argument).evalf(50).has(HPL)


def test_float_constructor_remains_lossless_symbolic():
    assert isinstance(HPL((1, 0), s.Float('0.3', 50)), HPL)
