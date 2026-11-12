import mpmath as mp
import pytest
import sympy as s

from hypra.hpl_numeric import hpl_evalf
from hypra.nielsen import NielsenPolylog, hpl_to_nielsen, nielsen_to_hpl, reduce_nielsen
from hypra.serialization import SerializationError, decode_expr, encode_expr
from hypra.symbolic import HPL, reduce_known_hpl

z = s.Symbol('z')


def test_lossless_notation_and_partial_reverse_conversion():
    expr = 3*NielsenPolylog(2,3,z)+HPL((1,0),z)
    assert hpl_to_nielsen(nielsen_to_hpl(expr)) == expr
    assert NielsenPolylog(2,3,z).rewrite(HPL) == HPL((0,0,1,1,1),z)
    for word in [(),(0,),(1,),(1,0),(0,-1),(0,1,0)]:
        assert hpl_to_nielsen(HPL(word,z)) == HPL(word,z)


def test_classical_reduction_is_opt_in():
    atom = NielsenPolylog(2,1,z)
    assert isinstance(atom,NielsenPolylog)
    assert reduce_nielsen(atom) == s.polylog(3,z)
    assert atom.rewrite(s.polylog) == s.polylog(3,z)
    assert reduce_nielsen(NielsenPolylog(1,2,z)) == NielsenPolylog(1,2,z)


@pytest.mark.parametrize('n,p',[(1,1),(1,3),(2,2),(3,1)])
def test_derivative_matches_hpl_kernel(n,p):
    atom = NielsenPolylog(n,p,z)
    actual = reduce_known_hpl(nielsen_to_hpl(s.diff(atom,z)))
    expected = reduce_known_hpl(s.diff(atom.as_hpl(),z))
    assert s.simplify(actual-expected) == 0
    assert s.diff(NielsenPolylog(n,p,z*z),z) == atom.fdiff(3).subs(z,z*z)*2*z


@pytest.mark.parametrize('n,p',[(1,2),(2,2),(3,1)])
def test_integral_normalization_independent_of_hpl_mapping(n,p):
    with mp.workdps(60):
        x = mp.mpf(1)/3
        integral = (-1)**(n+p-1)/(mp.factorial(n-1)*mp.factorial(p))*mp.quad(
            lambda t: mp.log(t)**(n-1)*mp.log1p(-x*t)**p/t,[0,1])
        value = hpl_evalf((0,)*n+(1,)*p,s.Rational(1,3),dps=45)
        assert abs(value-integral) < mp.mpf('1e-45')


def test_ast_latex_and_invalid_indices():
    atom = NielsenPolylog(2,3,z)
    assert decode_expr(encode_expr(atom)) == atom
    assert s.latex(atom) == r'S_{2,3}\!\left(z\right)'
    for n,p in [(0,1),(1,0),(True,1),(2.0,3),(s.Symbol('n'),1)]:
        with pytest.raises(ValueError):
            NielsenPolylog(n,p,z)
    ast = encode_expr(atom)
    ast['args'][0] = {'head':'Integer','value':'0'}
    with pytest.raises(SerializationError):
        decode_expr(ast)
    with pytest.raises(TypeError):
        nielsen_to_hpl('NielsenPolylog(1,1,z)')
