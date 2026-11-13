import mpmath as mp
import pytest
import sympy as s

from hypra import NativeConfig, UnsupportedExpressionError, expand_epsilon, hypergeometric

e = s.Symbol("epsilon", real=True)
half = s.S.Half


def expand(upper, lower, through, assumptions=None):
    return expand_epsilon(
        hypergeometric(upper, lower, 1), e, through=through,
        assumptions=assumptions, config=NativeConfig(use_cache=False),
    )


def test_matched_half_rank_three():
    actual = expand([half, 1, 1+e], [half+1, 2], 1)
    # Integrate the binomial series for (1-x*t**2)**(-1-e) on the unit square.
    expected = 2*s.log(2) + e*(s.zeta(2)-2*s.log(2)**2)
    assert s.simplify(actual.polynomial-expected) == 0


def test_matched_half_rank_four_with_external_pole():
    h = hypergeometric([half, 1, 1, 1+e], [half+1, 2, 2], 1)
    actual = expand_epsilon(h/e, e, through=0, config=NativeConfig(use_cache=False))
    assert actual.coefficient(-1) == 4*s.log(2)-s.zeta(2)
    # The derivative of (1+e)_n/n! is H_n. Independently sum the
    # resulting series with a bound on the omitted positive terms.
    with mp.workdps(30):
        harmonic = mp.mpf(0)
        expected = mp.mpf(0)
        for n in range(1, 2001):
            harmonic += mp.mpf(1)/n
            expected += harmonic/((2*n+1)*(n+1)**2)
        value = mp.mpf(str(actual.coefficient(0).evalf(20)))
        assert 0 < value-expected < mp.mpf("1e-6")


def test_unmatched_upper_half_rank_three():
    actual = expand([half, e, e], [1, 1], 2)
    assert s.simplify(actual.polynomial-1-e**2*(s.zeta(2)-2*s.log(2)**2)) == 0


def test_unmatched_lower_half_rank_three():
    actual = expand([e, e, e], [half, 1], 3)
    # The first nontrivial coefficient is sum 4**n/(n**3 binomial(2*n,n)).
    with mp.workdps(30):
        expected = mp.quad(lambda t: 2*mp.polylog(2, 1-t*t)/(1-t*t), [0, 1])
        value = mp.mpf(str(actual.coefficient(3).evalf(22)))
        assert abs(value-expected) < mp.mpf("1e-19")


def test_critical_rank_three_requires_convergent_side_and_plans_pole():
    upper, lower = [1, 1, 1], [2, 1+e]
    # The Euler beta integral gives epsilon*polygamma(1, epsilon).
    result = expand(upper, lower, 1, e > 0)
    assert s.simplify(result.polynomial-1/e-e*s.zeta(2)) == 0
    for assumptions in (None, e < 0):
        with pytest.raises(UnsupportedExpressionError, match="positive excess"):
            expand(upper, lower, 1, assumptions)


def test_critical_half_integer_rank_three():
    # The residue and finite term also follow by integrating the Euler beta
    # representation of 2F1(1,1/2;1/2+epsilon;z) in z from zero to one.
    result = expand([half, 1, 1], [half+e, 2], 0, e > 0)
    assert s.simplify(result.polynomial - (1/e-2*s.log(2))) == 0
