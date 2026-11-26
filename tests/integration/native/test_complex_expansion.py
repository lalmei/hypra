import pytest
import sympy as s

from hypra import NativeConfig, UnsupportedExpressionError, expand_epsilon

e = s.Symbol("eps")


@pytest.mark.parametrize("z", [s.I / 3, (1 + s.I) / 4, (1 - s.I) / 4])
def test_explicit_complex_disk_exact_p_family(z):
    result = expand_epsilon(
        s.hyper([1, -e], [1 - 2 * e], z),
        e,
        2,
        config=NativeConfig(argument_domain="principal-disk", use_cache=False),
    )
    expected = 1 + e * s.log(1 - z) + e**2 * (-2 * s.polylog(2, z) - s.log(1 - z) ** 2 / 2)
    assert s.simplify(result.polynomial - expected) == 0
    assert result.provenance["backend"]["argument_domain"] == "principal-disk"


@pytest.mark.parametrize("z", [-1, 2 * s.I, s.Symbol("z")])
def test_disk_contract_rejects_unproved_or_outside_points(z):
    with pytest.raises(UnsupportedExpressionError, match="Principal-disk"):
        expand_epsilon(
            s.hyper([1, -e], [1 - 2 * e], z),
            e,
            1,
            config=NativeConfig(argument_domain="principal-disk"),
        )


def test_default_does_not_widen_complex_domain():
    with pytest.raises(UnsupportedExpressionError, match="0 <"):
        expand_epsilon(s.hyper([1, -e], [1 - 2 * e], s.I / 3), e, 1)


@pytest.mark.parametrize("z", [3 * s.I, 2 + s.I, 2 - s.I])
def test_explicit_cut_plane_closed_forms(z):
    result = expand_epsilon(
        s.hyper([1, -e], [1 - 2 * e], z),
        e,
        2,
        config=NativeConfig(argument_domain="principal-cut-plane", use_cache=False),
    )
    expected = 1 + e * s.log(1 - z) + e**2 * (-2 * s.polylog(2, z) - s.log(1 - z) ** 2 / 2)
    assert s.simplify(result.polynomial - expected) == 0
    # Secondary numerical check at a small nonzero regulator, not the derivation.
    numeric = s.N(result.polynomial.subs(e, s.Rational(1, 100000)), 35)
    oracle = s.N(s.hyper([1, -s.Rational(1, 100000)], [1 - s.Rational(2, 100000)], z), 35)
    assert abs(complex(numeric - oracle)) < 1e-12


@pytest.mark.parametrize("z", [1, 2, s.Symbol("z")])
def test_cut_plane_requires_proved_nonboundary_point(z):
    with pytest.raises(UnsupportedExpressionError):
        expand_epsilon(
            s.hyper([1, -e], [1 - 2 * e], z),
            e,
            1,
            config=NativeConfig(argument_domain="principal-cut-plane"),
        )


@pytest.mark.parametrize("z", [-s.S.Half, -1, -3, -s.sqrt(2)])
@pytest.mark.parametrize("domain", ["real-interval", "principal-cut-plane"])
def test_negative_real_integer_expansion(z, domain):
    result = expand_epsilon(
        s.hyper([1, e], [1 + e], z),
        e,
        2,
        config=NativeConfig(argument_domain=domain, use_cache=False),
    )
    # (e)_n/(1+e)_n=e/(n+e). The germ continues along the entire negative ray.
    assert s.simplify(result.polynomial - (1 - e * s.log(1 - z) - e**2 * s.polylog(2, z))) == 0


def test_principal_disk_includes_negative_interior_point():
    result = expand_epsilon(
        s.hyper([1, e], [1 + e], -s.S.Half),
        e,
        1,
        config=NativeConfig(argument_domain="principal-disk", use_cache=False),
    )
    assert s.simplify(result.coefficient(1) + s.log(s.Rational(3, 2))) == 0


def test_negative_symbolic_domain_preserves_exact_mixed_integer_hpl():
    import mpmath as mp

    from hypra import HPL

    z = s.Symbol("z")
    result = expand_epsilon(
        s.hyper([1, e], [1 - 2 * e], z), e, 3, z < 0, config=NativeConfig(use_cache=False)
    )
    # The integer word (0,1,1) survives classical reduction at this order.
    mixed = result.coefficient(3).atoms(HPL)
    assert mixed and all(s.cancel(atom.argument - z / (z - 1)) == 0 for atom in mixed)
    with mp.workdps(45):
        expected = mp.diff(lambda t: mp.hyp2f1(1, t, 1 - 2 * t, -2), 0, 3) / 6
        actual = result.coefficient(3).subs(z, -2).evalf(32)
        assert abs(actual - s.Float(str(expected), 40)) < s.Rational(1, 10**27)


def test_negative_moving_lower_pole_with_symbolic_slope():
    a = s.Symbol("a")
    f = s.hyper([1, 1 + a * e], [a * e], -2)
    result = expand_epsilon(f, e, 0, s.Ne(a, 0), config=NativeConfig(use_cache=False))
    # (1+a*e)_n/(a*e)_n = 1 + n/(a*e). Sum the geometric series and its derivative.
    assert s.simplify(result.polynomial - (s.S.One / 3 - 2 / (9 * a * e))) == 0


@pytest.mark.parametrize(
    "upper,lower,order",
    [
        ([e, 2 * e], [s.S.Half + 3 * e], 3),
        ([s.S.Half + e, 2 * e], [1 + 3 * e], 3),
        ([1 + e, s.S.Half + e], [s.S.Half + 2 * e], 2),
        (
            [s.S.Half + e, s.S.Half + 2 * e, s.S.Half + 3 * e],
            [s.S.Half + 4 * e, s.S.Half + 5 * e],
            2,
        ),
    ],
    ids=["single-half-lower", "single-half-upper", "matched", "all-half"],
)
def test_negative_real_half_families_against_defining_function(upper, lower, order):
    import mpmath as mp

    f = s.hyper(upper, lower, -2, evaluate=False)
    result = expand_epsilon(f, e, order, config=NativeConfig(use_cache=False))
    # Differentiate the independent arbitrary-precision hypergeometric implementation.
    evaluate = s.lambdify(e, f, "mpmath")
    with mp.workdps(45):
        for k in range(order + 1):
            expected = mp.diff(evaluate, 0, k) / mp.factorial(k)
            difference = result.coefficient(k).evalf(32) - s.Float(str(expected), 40)
            assert abs(difference) < s.Rational(1, 10**27)
