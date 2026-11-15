import mpmath as mp
import pytest
import sympy as s

from hypra import NativeConfig, UnsupportedExpressionError, expand_epsilon, hypergeometric

e = s.Symbol("epsilon", real=True)
half = s.S.Half


def expand(upper, lower, through, assumptions=None):
    return expand_epsilon(
        hypergeometric(upper, lower, 1),
        e,
        through=through,
        assumptions=assumptions,
        config=NativeConfig(use_cache=False),
    )


def test_matched_half_rank_three():
    actual = expand([half, 1, 1 + e], [half + 1, 2], 1)
    # Integrate the binomial series for (1-x*t**2)**(-1-e) on the unit square.
    expected = 2 * s.log(2) + e * (s.zeta(2) - 2 * s.log(2) ** 2)
    assert s.simplify(actual.polynomial - expected) == 0


def test_matched_half_rank_four_with_external_pole():
    h = hypergeometric([half, 1, 1, 1 + e], [half + 1, 2, 2], 1)
    actual = expand_epsilon(h / e, e, through=0, config=NativeConfig(use_cache=False))
    assert actual.coefficient(-1) == 4 * s.log(2) - s.zeta(2)
    # The derivative of (1+e)_n/n! is H_n. Independently sum the
    # resulting series with a bound on the omitted positive terms.
    with mp.workdps(30):
        harmonic = mp.mpf(0)
        expected = mp.mpf(0)
        for n in range(1, 2001):
            harmonic += mp.mpf(1) / n
            expected += harmonic / ((2 * n + 1) * (n + 1) ** 2)
        value = mp.mpf(str(actual.coefficient(0).evalf(20)))
        assert 0 < value - expected < mp.mpf("1e-6")


def test_unmatched_upper_half_rank_three():
    actual = expand([half, e, e], [1, 1], 2)
    assert s.simplify(actual.polynomial - 1 - e**2 * (s.zeta(2) - 2 * s.log(2) ** 2)) == 0


def test_unmatched_lower_half_rank_three():
    actual = expand([e, e, e], [half, 1], 3)
    # The first nontrivial coefficient is sum 4**n/(n**3 binomial(2*n,n)).
    with mp.workdps(30):
        expected = mp.quad(lambda t: 2 * mp.polylog(2, 1 - t * t) / (1 - t * t), [0, 1])
        value = mp.mpf(str(actual.coefficient(3).evalf(22)))
        assert abs(value - expected) < mp.mpf("1e-19")


def test_critical_rank_three_requires_convergent_side_and_plans_pole():
    upper, lower = [1, 1, 1], [2, 1 + e]
    # The Euler beta integral gives epsilon*polygamma(1, epsilon).
    result = expand(upper, lower, 1, e > 0)
    assert s.simplify(result.polynomial - 1 / e - e * s.zeta(2)) == 0
    for assumptions in (None, e < 0):
        with pytest.raises(UnsupportedExpressionError, match="positive excess"):
            expand(upper, lower, 1, assumptions)


def test_critical_half_integer_rank_three():
    # The residue and finite term also follow by integrating the Euler beta
    # representation of 2F1(1,1/2;1/2+epsilon;z) in z from zero to one.
    result = expand([half, 1, 1], [half + e, 2], 0, e > 0)
    assert s.simplify(result.polynomial - (1 / e - 2 * s.log(2))) == 0


def test_meromorphic_negative_integer_excess_with_moving_lower_pole():
    h = hypergeometric([1, 1, 1], [2, e], 1)
    config = NativeConfig(unit_endpoint="meromorphic", use_cache=False)
    result = expand_epsilon(h, e, through=-2, config=config)
    # On Re(e)>1 the Euler beta integral gives (e-1)*trigamma(e-1).
    assert result.polynomial == -1 / e**2
    assert result.provenance["unit_endpoint"] == "meromorphic"
    with pytest.raises(UnsupportedExpressionError, match="positive excess"):
        expand_epsilon(h, e, through=-2)


def test_meromorphic_half_integer_rank_three_reduces_by_factorial_moments():
    # Removing (3/2)_n/(1/2)_n = 2*n+1 leaves Gauss sums.
    h = hypergeometric([1, 2 + e, s.Rational(3, 2)], [2, half], 1)
    result = expand_epsilon(h, e, through=0, config=NativeConfig(unit_endpoint="meromorphic"))
    # The two Gauss quotients combine to 1/(1+epsilon).
    assert result.polynomial == 1


def test_meromorphic_half_integer_rank_four():
    h = hypergeometric([1, 2, 3, half], [2, 2, s.Rational(3, 2) + e], 1)
    result = expand_epsilon(h, e, through=0, config=NativeConfig(unit_endpoint="meromorphic"))
    # Cancel 2, use (3)_n/(2)_n=(n+2)/2, then apply Gauss to each moment.
    a, b, c = s.S.One, half, s.Rational(3, 2) + e
    gauss = lambda a, b, c: s.gamma(c) * s.gamma(c - a - b) / (s.gamma(c - a) * s.gamma(c - b))
    reference = gauss(a, b, c) + a * b / (2 * c) * gauss(a + 1, b + 1, c + 1)
    assert s.simplify(result.polynomial - s.series(reference, e, 0, 1).removeO()) == 0


def test_meromorphic_fixed_excess_pole_rejected():
    h = hypergeometric([1, 1, 1], [half, half], 1)
    with pytest.raises(UnsupportedExpressionError, match="fixed excess pole"):
        expand_epsilon(h, e, through=0, config=NativeConfig(unit_endpoint="meromorphic"))


def test_excess_recurrence_agrees_with_independent_convergent_series():
    from hypra.native.endpoints import rewrite_endpoints

    h = hypergeometric([half, 1, 1], [half - 1 + e, 2], 1)
    rewritten = rewrite_endpoints(h, e, meromorphic=True)
    assert all(sum(atom.bq).subs(e, 0) > sum(atom.ap) for atom in rewritten.atoms(s.hyper))
    point = s.Rational(5, 2)
    assert abs((rewritten.subs(e, point) - h.subs(e, point)).evalf(22)) < s.Float("1e-20")


def test_fixed_gauss_excess_transverse_laurent_series():
    config = NativeConfig(unit_endpoint="meromorphic", unit_excess_shift=1)
    h = hypergeometric([1, 1], [2], 1)
    result = expand_epsilon(h, e, through=0, config=config)
    assert result.polynomial == 1 / e + 1
    assert result.provenance["unit_excess_shift"] == 1


def test_fixed_negative_half_excess_transverse_rank_three():
    h = hypergeometric([half, 1, 1], [-half, 2], 1)
    config = NativeConfig(unit_endpoint="meromorphic", unit_excess_shift=2)
    result = expand_epsilon(h, e, through=0, config=config)
    # (1/2)_n/(-1/2)_n = 1-2*n. Its two factorial moments reduce
    # to (1+e)/e - 2*(1+e)/(e*(e-1)).
    reference = (1 + e) / e - 2 * (1 + e) / (e * (e - 1))
    assert result.polynomial == s.series(reference, e, 0, 1).removeO()


def test_fixed_excess_prescription_bounds_and_validation():
    h = hypergeometric([1, 1], [2], 1)
    with pytest.raises(ValueError, match="meromorphic mode"):
        NativeConfig(unit_excess_shift=1)
    with pytest.raises(UnsupportedExpressionError, match="reduced lower parameter count"):
        expand_epsilon(
            h, e, through=0, config=NativeConfig(unit_endpoint="meromorphic", unit_excess_shift=2)
        )


def test_fixed_negative_half_excess_transverse_rank_four():
    h = hypergeometric([half, 1, 1, 3], [-half, 2, 2], 1)
    config = NativeConfig(unit_endpoint="meromorphic", unit_excess_shift=3)
    result = expand_epsilon(h, e, through=0, config=config)
    # The unshifted ratios multiply to (1-2*n)*(n+2)/2.
    # Sum the three factorial moments of the remaining Gauss series.
    reference = (1 + e) / e * (1 - s.Rational(5, 2) / (e - 1) - 4 / ((e - 1) * (e - 2)))
    assert result.polynomial == s.series(reference, e, 0, 1).removeO()
