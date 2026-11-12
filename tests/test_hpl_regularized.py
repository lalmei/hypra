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
