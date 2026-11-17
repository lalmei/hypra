"""Named-function reductions checked by derivatives, origins and integrals."""

import mpmath as mp
import pytest
import sympy as s

from hypra import HPL, NielsenPolylog, UnsupportedExpressionError, transform_polylog_arguments
from hypra.hpl.algebra import shuffle_expand
from hypra.polylog_reduction import reduce_polylogs

Z = s.Symbol("z", positive=True)
ARGUMENTS = [Z, 1 - Z, Z / (Z - 1), 1 - 1 / Z, 1 / Z, 1 / (1 - Z)]


@pytest.mark.parametrize("n,p", [(2, 2), (2, 3)])
def test_same_argument_named_reduction_has_exact_derivative_and_origin(n, p):
    source = NielsenPolylog(n, p, Z)
    result = reduce_polylogs(source, Z, domain="0<z<1")
    assert not result.has(HPL, source)
    assert result.atoms(NielsenPolylog) <= {NielsenPolylog(3, 2, Z)}
    derivative = transform_polylog_arguments(s.diff(result - source, Z), Z, domain="0<z<1")
    assert s.simplify(shuffle_expand(derivative, Z)) == 0
    # From its defining integral, S32(z)=z**2/16+O(z**3). Replacing this
    # bounded, vanishing summand by its first term preserves the origin limit
    # and avoids SymPy's generic Function limit attempting 0/0 derivatives.
    origin_germ = result.xreplace({NielsenPolylog(3, 2, Z): Z**2 / 16})
    assert s.limit(origin_germ, Z, 0, dir="+") == 0
    assert reduce_polylogs(result, Z, domain="0<z<1") == result


def _defining_integral(n, p, argument, sign):
    def integrand(t):
        residual = 1 - argument * t
        if not t or not residual:
            return mp.mpf(0)
        logarithm = mp.log(residual) if residual > 0 else mp.log(-residual) - sign * mp.pi * 1j
        return mp.log(t) ** (n - 1) * logarithm**p / t

    points = [0, 1 / argument, 1] if argument > 1 else [0, 1]
    return (
        (-1) ** (n + p - 1) * mp.quad(integrand, points) / (mp.factorial(n - 1) * mp.factorial(p))
    )


@pytest.mark.parametrize("n,p", [(2, 2), (2, 3), (3, 2)])
@pytest.mark.parametrize("argument", ARGUMENTS)
@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
def test_all_nielsen_cross_ratios_in_named_basis_against_integrals(n, p, argument, boundary, sign):
    result = reduce_polylogs(NielsenPolylog(n, p, argument), Z, domain="0<z<1", boundary=boundary)
    assert not result.has(HPL)
    assert result.atoms(NielsenPolylog) <= {NielsenPolylog(3, 2, Z)}
    with mp.workdps(45):
        x = mp.mpf(str(argument.subs(Z, s.Rational(2, 5)).evalf(48)))
        expected = _defining_integral(n, p, x, sign)
        value = result.subs(Z, s.Rational(2, 5)).evalf(40)
        computed = mp.mpc(str(s.re(value)), str(s.im(value)))
        assert abs(computed - expected) < mp.mpf("1e-34")


@pytest.mark.parametrize("n,p", [(1, 2), (1, 3), (1, 4)])
def test_remaining_low_weight_nielsen_words_share_the_basis(n, p):
    result = reduce_polylogs(NielsenPolylog(n, p, 1 - Z), Z, domain="0<z<1")
    assert not result.has(HPL)
    assert result.atoms(NielsenPolylog) <= {NielsenPolylog(3, 2, Z)}
    with mp.workdps(45):
        expected = _defining_integral(n, p, mp.mpf(3) / 5, 1)
        assert abs(mp.mpf(str(result.subs(Z, s.Rational(2, 5)).evalf(40))) - expected) < mp.mpf(
            "1e-34"
        )


@pytest.mark.parametrize("weight", [2, 3, 4, 5])
@pytest.mark.parametrize("boundary,sign", [("upper", 1), ("lower", -1)])
def test_classical_inversions_keep_the_requested_cut_side(weight, boundary, sign):
    result = reduce_polylogs(s.polylog(weight, 1 / Z), Z, domain="0<z<1", boundary=boundary)
    assert not result.has(HPL, NielsenPolylog)
    with mp.workdps(45):
        expected = _defining_integral(weight - 1, 1, mp.mpf(5) / 2, sign)
        value = result.subs(Z, s.Rational(2, 5)).evalf(40)
        computed = mp.mpc(str(s.re(value)), str(s.im(value)))
        assert abs(computed - expected) < mp.mpf("1e-34")


def test_canonical_s32_and_endpoint_constant():
    canonical = NielsenPolylog(3, 2, Z)
    assert reduce_polylogs(canonical, Z, domain="0<z<1") == canonical
    actual = reduce_polylogs(NielsenPolylog(3, 2, -1), Z, domain="0<z<1")
    assert s.simplify(actual - s.pi**2 * s.zeta(3) / 12 + 29 * s.zeta(5) / 32) == 0


@pytest.mark.parametrize("argument", ARGUMENTS[-2:])
def test_positive_cut_requires_explicit_side(argument):
    with pytest.raises(UnsupportedExpressionError, match="boundary"):
        reduce_polylogs(NielsenPolylog(2, 3, argument), Z, domain="0<z<1")


def test_unrecognized_families_and_unrelated_factors_are_retained():
    untouched = s.polylog(7, 1 / Z) + NielsenPolylog(2, 4, Z) + NielsenPolylog(2, 2, 2 * Z)
    assert reduce_polylogs(untouched, Z, domain="0<z<1") == untouched
    scalar = s.Symbol("a")
    atom = NielsenPolylog(2, 2, Z)
    assert reduce_polylogs(scalar / Z * atom, Z, domain="0<z<1") == (
        scalar / Z * reduce_polylogs(atom, Z, domain="0<z<1")
    )


def test_domain_exactness_and_argument_contracts():
    atom = NielsenPolylog(2, 2, Z)
    with pytest.raises(UnsupportedExpressionError, match="domain"):
        reduce_polylogs(atom, Z, domain="complex")
    with pytest.raises(ValueError, match="boundary"):
        reduce_polylogs(atom, Z, domain="0<z<1", boundary="unspecified")
    with pytest.raises(UnsupportedExpressionError, match="exact finite"):
        reduce_polylogs(atom + s.Float("0.1"), Z, domain="0<z<1")
    with pytest.raises(TypeError, match="symbolic"):
        reduce_polylogs(atom, Z + 1, domain="0<z<1")
    with pytest.raises(TypeError, match="commutative"):
        reduce_polylogs(atom, s.Symbol("a", commutative=False), domain="0<z<1")
