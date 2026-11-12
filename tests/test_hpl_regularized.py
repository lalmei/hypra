import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import shuffle_words
from hypra.hpl_alternating import hpl_alternating_at_one
from hypra.hpl_regularized import hpl_endpoint_polynomial as polynomial
from hypra.hpl_regularized import shuffle_regularized_hpl_at_one as regularized

T = s.Symbol('T')


def test_explicit_divergent_polynomials_not_finite_limits():
    assert polynomial((1,), T) == T
    assert regularized((1,)) == 0
    assert polynomial((1,1,1), T) == T**3/6
    with pytest.raises(UnsupportedExpressionError):
        hpl_alternating_at_one((1,))
    assert polynomial((1,0,1), T) == T*s.zeta(2)-2*s.zeta(3)
    assert regularized((1,0,1)) == -2*s.zeta(3)


@pytest.mark.parametrize('word', [(),(0,),(-1,),(0,1),(1,0),(1,0,0),(-1,0,1)])
def test_convergent_words_have_existing_endpoint(word):
    assert s.simplify(polynomial(word,T)-hpl_alternating_at_one(word)) == 0


@pytest.mark.parametrize('left,right', [((1,),(-1,)),((1,),(1,0)),((1,1),(0,1))])
def test_shuffle_product_regularization(left,right):
    shuffled = sum(count*polynomial(word,T) for word,count in shuffle_words(left,right).items())
    assert s.expand(shuffled-polynomial(left,T)*polynomial(right,T)) == 0
    constant = sum(count*regularized(word) for word,count in shuffle_words(left,right).items())
    assert s.expand(constant-regularized(left)*regularized(right)) == 0


def test_explicit_symbol_and_word_validation():
    for bad in ['T',1,s.Symbol('x')+1]:
        with pytest.raises(TypeError):
            polynomial((1,),bad)
    with pytest.raises(UnsupportedExpressionError):
        regularized((2,))


def test_regularization_resource_budget():
    with pytest.raises(UnsupportedExpressionError, match="weight-12"):
        polynomial((1,)*13, s.Symbol("T"))


def test_expression_endpoint_cancellations_and_products():
    from hypra.hpl_regularized import hpl_expression_endpoint_polynomial as endpoint
    from hypra.symbolic import HPL

    z, a = s.symbols('z a')
    h1, hm = HPL((1,),z), HPL((-1,),z)
    assert endpoint(h1**2-2*HPL((1,1),z),z,T) == 0
    assert endpoint(a*h1**2+hm,z,T) == a*T**2+s.log(2)
    expr = h1*hm-HPL((1,-1),z)
    assert s.expand(endpoint(expr,z,T)-hpl_alternating_at_one((-1,1))) == 0
    assert endpoint(s.Integer(3),z,T) == 3


@pytest.mark.parametrize('kind', ['inverse','pole','regular_z','function','argument','float','T','power'])
def test_expression_endpoint_outside_ring_rejected(kind):
    from hypra.hpl_regularized import hpl_expression_endpoint_polynomial as endpoint
    from hypra.symbolic import HPL

    z = s.Symbol('z')
    h = HPL((1,),z)
    expressions = {'inverse':1/h, 'pole':h/(1-z), 'regular_z':z*h,
                   'function':s.sin(h), 'argument':HPL((1,),1-z),
                   'float':s.Float('0.5')*h, 'T':T+h, 'power':h**13}
    with pytest.raises(UnsupportedExpressionError):
        endpoint(expressions[kind],z,T)


def test_expression_endpoint_validates_symbols_and_text():
    from hypra.hpl_regularized import hpl_expression_endpoint_polynomial as endpoint
    z = s.Symbol('z')
    with pytest.raises(TypeError):
        endpoint('HPL((1,),z)',z,T)
    with pytest.raises(TypeError):
        endpoint(s.Integer(1),z,z)
