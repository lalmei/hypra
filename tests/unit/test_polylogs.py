"""Independent defining-integral checks of the six polylog argument maps."""

import mpmath as mp
import pytest
import sympy as s

from hypra import HPL, NielsenPolylog, UnsupportedExpressionError, transform_polylog_arguments
from hypra.hpl.algebra import shuffle_expand

Z = s.Symbol("z", positive=True)
ARGUMENTS = [Z, 1 - Z, Z / (Z - 1), 1 - 1 / Z, 1 / Z, 1 / (1 - Z)]
INDICES = [(1, 1), (2, 1), (3, 1), (4, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 2)]


def _defining_integral(n, p, argument, sign):
    # The path is the positive integration interval with the requested exact
    # boundary value of log(1-x*t); no artificial imaginary offset is used.
    def integrand(t):
        if not t:
            return mp.mpf(0)
        residual = 1 - argument * t
        if not residual:
            return mp.mpf(0)  # Single point; quadrature uses open subintervals.
        logarithm = mp.log(residual) if residual > 0 else mp.log(-residual) - sign * mp.pi * 1j
        return mp.log(t) ** (n - 1) * logarithm**p / t

    points = [0, 1 / argument, 1] if argument > 1 else [0, 1]
    return (
        (-1) ** (n + p - 1) * mp.quad(integrand, points) / (mp.factorial(n - 1) * mp.factorial(p))
    )


@pytest.mark.parametrize("argument", ARGUMENTS)
@pytest.mark.parametrize("n,p", INDICES)
@pytest.mark.parametrize("side,sign", [("upper", 1), ("lower", -1)])
def test_appendix_argument_families_against_defining_integrals(argument, n, p, side, sign):
    source = s.polylog(n + 1, argument) if p == 1 else NielsenPolylog(n, p, argument)
    result = transform_polylog_arguments(source, Z, domain="0<z<1", boundary=side)
    assert all(atom.argument == Z for atom in result.atoms(HPL))
    with mp.workdps(45):
        expected = _defining_integral(
            n, p, mp.mpf(str(argument.subs(Z, s.Rational(2, 5)).evalf(48))), sign
        )
        actual = result.subs(Z, s.Rational(2, 5)).evalf(40)
        real, imag = actual.as_real_imag()
        value = mp.mpc(str(real), str(imag))
        assert abs(value - expected) < mp.mpf("1e-34")


def test_landen_identity_exactly():
    result = transform_polylog_arguments(s.polylog(2, Z / (Z - 1)), Z, domain="0<z<1")
    assert s.expand(result + s.polylog(2, Z) + s.log(1 - Z) ** 2 / 2) == 0


def test_reflection_identity_exactly():
    source = s.polylog(2, Z) + s.polylog(2, 1 - Z)
    result = transform_polylog_arguments(source, Z, domain="0<z<1")
    assert shuffle_expand(result + s.log(Z) * s.log(1 - Z) - s.pi**2 / 6, Z) == 0


@pytest.mark.parametrize("order", [2, 3, 4, 5, 8])
def test_classical_inversion_jump(order):
    source = s.polylog(order, 1 / Z)
    upper = transform_polylog_arguments(source, Z, domain="0<z<1", boundary="upper")
    lower = transform_polylog_arguments(source, Z, domain="0<z<1", boundary="lower")
    jump = 2 * s.pi * s.I * (-s.log(Z)) ** (order - 1) / s.factorial(order - 1)
    assert s.simplify(upper - lower - jump) == 0


@pytest.mark.parametrize("argument", ARGUMENTS[-2:])
def test_positive_cut_requires_side(argument):
    with pytest.raises(UnsupportedExpressionError, match="boundary"):
        transform_polylog_arguments(s.polylog(2, argument), Z, domain="0<z<1")


def test_constant_arguments_and_unrelated_terms_are_preserved():
    expression = Z**2 + NielsenPolylog(3, 2, -1) + s.polylog(3, s.S.Half)
    assert transform_polylog_arguments(expression, Z, domain="0<z<1") == expression


@pytest.mark.parametrize(
    "expression", [s.polylog(2, 2 * Z), s.polylog(9, 1 - Z), s.polylog(s.Rational(3, 2), Z)]
)
def test_unrecognized_arguments_and_orders_rejected(expression):
    with pytest.raises(UnsupportedExpressionError):
        transform_polylog_arguments(expression, Z, domain="0<z<1")
