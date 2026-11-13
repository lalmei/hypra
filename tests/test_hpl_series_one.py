import itertools

import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_numeric import hpl_evalf
from hypra.hpl_series_one import hpl_series_at_one

t = s.Symbol('t', positive=True)


@pytest.mark.parametrize('word,exact', [
    ((),s.S.One), ((0,),s.log(1-t)), ((1,),-s.log(t)), ((-1,),s.log(2-t)),
    ((0,0),s.log(1-t)**2/2), ((1,1),s.log(t)**2/2),
    ((-1,-1),s.log(2-t)**2/2),
])
def test_elementary_exact_series(word,exact):
    actual = hpl_series_at_one(word,t,through=4).removeO()
    expected = exact.series(t,0,5).removeO()
    assert s.expand(actual-expected) == 0


def test_dilog_near_one_independent_reflection():
    expected = s.pi**2/6-s.log(t)*s.log(1-t)-s.polylog(2,t)
    actual = hpl_series_at_one((0,1),t,through=4)
    assert s.expand(actual.removeO()-expected.series(t,0,5).removeO()) == 0
    assert actual.getO() == s.Order(t**5*s.log(t)**2,t)


def test_all_ternary_weight_two_differential_recurrence():
    for word in itertools.product((-1,0,1),repeat=2):
        actual = hpl_series_at_one(word,t,through=3).removeO()
        tail = hpl_series_at_one(word[1:],t,through=3).removeO()
        kernel = {1:-1/t,0:-1/(1-t),-1:-1/(2-t)}[word[0]]
        residual = s.diff(actual,t)-kernel*tail
        assert s.expand(residual.series(t,0,3).removeO()) == 0


@pytest.mark.parametrize('word', [(1,0,1),(-1,0,1),(0,-1,0)])
def test_mixed_word_numerical_secondary(word):
    with mp.workdps(50):
        delta = s.Rational(1,50)
        expansion = hpl_series_at_one(word,t,through=5).removeO().subs(t,delta)
        approximation = mp.mpf(str(expansion.evalf(45)))
        reference = hpl_evalf(word,1-delta,dps=40)
        assert abs(approximation-reference) < mp.mpf('1e-8')


def test_limits_and_no_false_plain_power_remainder():
    value = hpl_series_at_one((0,1),t,through=0)
    assert value.getO() == s.Order(t*s.log(t)**2,t)
    for bad in [-1,25,True]:
        with pytest.raises(UnsupportedExpressionError):
            hpl_series_at_one((0,1),t,through=bad)
    with pytest.raises(UnsupportedExpressionError):
        hpl_series_at_one((1,)*9,t)
    with pytest.raises(TypeError):
        hpl_series_at_one((0,1),'t')
