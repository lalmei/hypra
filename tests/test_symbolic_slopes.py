import pytest
import sympy as s

from hypra import HPL, expand_epsilon
from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import hpl_series
from hypra.planning import _affine_rational

eps, z, a, b, c = s.symbols('eps z a b c')
assumptions = s.And(z > 0, z < 1)


@pytest.mark.parametrize('upper,lower', [
    ([1+a*eps,b*eps],[1+c*eps]),
    ([2+a*eps,b*eps],[3+c*eps]),
    ([s.Rational(1,2)+a*eps,b*eps],[s.Rational(1,2)+c*eps]),
])
def test_symbolic_coefficients_and_degenerate_specializations(upper,lower):
    h = s.hyper(upper,lower,z,evaluate=False)
    result = expand_epsilon(h,eps,through=2,assumptions=assumptions).polynomial
    assert result.is_polynomial(a,b,c)
    for substitutions in [{a:0,b:0,c:0},{a:1,b:-1,c:1},{a:2,b:1,c:0}]:
        direct = expand_epsilon(h.subs(substitutions),eps,through=2,assumptions=assumptions).polynomial
        difference = s.simplify(result.subs(substitutions)-direct)
        if difference != 0:
            # Independent exact local coefficients avoid assuming log-product
            # identities across branches or relying on a complete HPL reducer.
            t = s.Symbol('t', positive=True)
            local = difference.subs(z,t**2)
            local = local.xreplace({h: hpl_series(h,t,8) for h in local.atoms(HPL)})
            assert s.series(local,t,0,8).removeO().expand() == 0


def test_polynomial_slope_and_exact_low_order_formula():
    result = expand_epsilon(s.hyper([1,(a*b+1)*eps],[1+c*eps],z),eps,
                            through=1,assumptions=assumptions).polynomial
    assert s.expand(result - (1-(a*b+1)*eps*s.log(1-z))) == 0


@pytest.mark.parametrize('parameter',[1+eps/a,1+eps*s.sin(a),1+eps*s.Float('0.5'),1+eps**2])
def test_nonpolynomial_inexact_or_nonlinear_slopes_rejected(parameter):
    with pytest.raises(UnsupportedExpressionError):
        _affine_rational(parameter,eps)


def test_symbolic_poles_and_argument_dependent_slopes_fail_closed():
    for h in [s.hyper([1+a*eps,b*eps],[c*eps],z,evaluate=False),
              s.hyper([1+z*eps,b*eps],[1+c*eps],z,evaluate=False)]:
        with pytest.raises(UnsupportedExpressionError):
            expand_epsilon(h,eps,through=1,assumptions=assumptions)


def test_symbolic_contiguous_removable_equal_slope_specialization():
    from hypra.contiguous import expand_contiguous
    from hypra.nested_sums import canonical_hyper_series

    def builder(bases, slopes, argument, regulator, depth):
        return canonical_hyper_series(list(zip(bases[0], slopes[0], strict=True)),
                                      slopes[1], argument, regulator, depth)
    result = expand_contiguous([(0,a),(1,b)],[(2,c)],[1,1],[1],z,eps,1,builder)
    assert result.is_polynomial(a,b,c)
    assert s.simplify(result.subs({a:0,b:0,c:0})-1) == 0


def test_moving_symbolic_lower_pole_requires_assumption_and_only_source_denominator():
    h = s.hyper([1,1+a*eps],[b*eps],z,evaluate=False)
    with pytest.raises(UnsupportedExpressionError, match='nonzero'):
        expand_epsilon(h,eps,through=1,assumptions=assumptions)
    approved = s.And(assumptions,s.Ne(b,0))
    result = expand_epsilon(h,eps,through=1,assumptions=approved).polynomial
    from hypra.contiguous import _hpl_pieces
    assert all(s.cancel(b*v).is_polynomial(a,b) for v in _hpl_pieces(result).values())
    for values in [{a:0,b:1},{a:1,b:1},{a:-1,b:2}]:
        direct = expand_epsilon(h.subs(values),eps,through=1,assumptions=assumptions).polynomial
        assert s.simplify(result.subs(values)-direct) == 0


def test_standalone_contiguous_lower_pole_requires_explicit_nonzero():
    from hypra.contiguous import expand_contiguous
    from hypra.nested_sums import canonical_hyper_series

    def builder(bases,slopes,argument,regulator,depth):
        return canonical_hyper_series(list(zip(bases[0],slopes[0],strict=True)),slopes[1],
                                      argument,regulator,depth)
    arguments = ([(1,0),(1,a)],[(0,b)],[1,1],[1],z,eps,1,builder)
    with pytest.raises(UnsupportedExpressionError, match='nonzero'):
        expand_contiguous(*arguments)
    result = expand_contiguous(*arguments,assumptions=s.Ne(b,0))
    from hypra.contiguous import _hpl_pieces
    assert all(s.cancel(b*v).is_polynomial(a,b) for v in _hpl_pieces(result).values())
