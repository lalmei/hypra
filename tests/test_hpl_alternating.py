import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_alternating import hpl_alternating_at_one as endpoint
from hypra.hpl_alternating import reduce_euler_sum
from hypra.hpl_endpoints import hpl_at_one
from hypra.serialization import SerializationError, decode_expr, encode_expr
from hypra.symbolic import EulerSum


def test_depth_one_and_repeated_letters():
    assert endpoint(()) == 1
    assert endpoint((0, 0)) == 0
    assert endpoint((-1,)) == s.log(2)
    assert endpoint((0, -1)) == s.zeta(2) / 2
    assert endpoint((0, 0, -1)) == 3 * s.zeta(3) / 4
    for n in range(1, 6):
        assert endpoint((-1,) * n) == s.log(2) ** n / s.factorial(n)
    assert endpoint((-1, 0)) == -s.zeta(2) / 2
    assert endpoint((1, 0)) == -s.zeta(2)


def test_mapping_and_trailing_zero_shuffle():
    assert endpoint((0, -1, 1), reduced=False) == -EulerSum((2, 1), (-1, -1))
    assert endpoint((-1, 0, -1), reduced=False) == EulerSum((1, 2), (-1, 1))
    assert endpoint((0, -1, 0), reduced=False) == 2 * EulerSum((3,), (-1,))
    assert endpoint((-1, 1, 0), reduced=False) == (
        EulerSum((2, 1), (-1, -1)) + EulerSum((1, 2), (-1, -1)))


@pytest.mark.parametrize('word', [(0, 1), (0, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1, 1)])
def test_non_alternating_matches_existing_endpoint(word):
    assert s.simplify(endpoint(word) - hpl_at_one(word)) == 0


@pytest.mark.parametrize('word', [(1,), (1, 1), (1, -1), (1, 0, -1), (1, -1, 0)])
def test_divergent_words_fail(word):
    with pytest.raises(UnsupportedExpressionError, match='Divergent'):
        endpoint(word)


def test_abel_endpoint_against_independent_nested_sum():
    # ES(2,1; -1,-1), computed directly with outer cutoff.
    nmax = 10000
    with mp.workdps(35):
        inner, value = mp.mpf(0), mp.mpf(0)
        for n in range(1, nmax + 1):
            value += (-1) ** n * inner / n**2
            inner += mp.mpf((-1) ** n) / n
        # H(0,-1,1;1) = integral H(-1,1;t)/t dt.
        # Swap finite positive kernel integrals to a single integral of logs.
        integral = mp.quad(lambda t: mp.log(t) * mp.log1p(-t) / (1 + t), [0, 1])
        assert abs(-value - integral) < mp.mpf('1e-8')


def test_lossless_ast_latex_and_no_float_cache_bypass():
    value = EulerSum((2, 1), (-1, -1))
    assert decode_expr(encode_expr(value)) == value
    assert 'operatorname{ES}' in s.latex(value)
    assert s.diff(value, s.Symbol('x')) == 0
    with pytest.raises(ValueError):
        EulerSum((2.0, 1), (-1, -1))
    with pytest.raises(ValueError):
        EulerSum((1,), (1,))
    with pytest.raises(ValueError):
        EulerSum((2,), ())
    ast = encode_expr(value)
    ast['args'][1]['args'][0] = {'head': 'Integer', 'value': '0'}
    with pytest.raises(SerializationError):
        decode_expr(ast)
    assert reduce_euler_sum(EulerSum((2, 2), (1, 1))) == 3 * s.zeta(4) / 4


def test_reduction_rejects_source_text():
    with pytest.raises(TypeError):
        reduce_euler_sum("1+1")


def test_repeated_color_newton_exact_low_depth():
    log2 = s.log(2)
    assert reduce_euler_sum(EulerSum((1, 1), (-1, -1))) == (log2**2 - s.zeta(2)) / 2
    expected = -log2**3 / 6 + log2 * s.zeta(2) / 2 - s.zeta(3) / 4
    assert s.expand(reduce_euler_sum(EulerSum((1,) * 3, (-1,) * 3)) - expected) == 0
    assert s.simplify(reduce_euler_sum(EulerSum((2,) * 3, (1,) * 3)) - s.pi**6 / 5040) == 0


@pytest.mark.parametrize('index,color,depth', [(1, -1, 4), (2, -1, 3), (2, 1, 4), (3, -1, 2)])
def test_newton_reduction_against_independent_finite_nested_sum(index, color, depth):
    with mp.workdps(40):
        # Descending update enumerates distinct indices directly, not Newton.
        coefficients = [mp.mpf(1)] + [mp.mpf(0)] * depth
        for n in range(1, 20001):
            term = mp.mpf(color**n) / n**index
            for j in range(depth, 0, -1):
                coefficients[j] += term * coefficients[j - 1]
        exact = reduce_euler_sum(EulerSum((index,) * depth, (color,) * depth))
        expected = mp.mpf(str(s.N(exact, 40)))
        assert abs(coefficients[depth] - expected) < mp.mpf('1e-4')
