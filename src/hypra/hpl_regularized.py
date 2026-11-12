"""Explicit shuffle-regularized real endpoint polynomials, not finite limits.

T=-log(1-z), z->1 from below. Polynomial constant terms define the shuffle
regularization. No stuffle regularization, path choice, or implicit replacement
of divergent endpoint values is performed.
"""

from functools import lru_cache

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import shuffle_words
from .hpl_alternating import hpl_alternating_at_one
from .hpl_numeric import _word


@lru_cache(maxsize=2048)
def _polynomial(word, T):
    if not word:
        return s.S.One
    if word[0] != 1:
        return hpl_alternating_at_one(word)
    if all(letter == 1 for letter in word):
        return T**len(word) / s.factorial(len(word))
    terms = shuffle_words((1,), word[1:])
    multiplicity = terms.pop(word)
    return s.expand((T * _polynomial(word[1:], T)
                     - sum(count * _polynomial(other, T) for other, count in terms.items()))
                    / multiplicity)


def hpl_endpoint_polynomial(word, T):
    """Return the polynomial asymptote in T=-log(1-z), z->1-, for one word.

    The omitted remainder tends to zero even after multiplication by a fixed
    power of T. Expanded letters are 0,+1,-1; weight limit is 12. For divergent
    words this is a regularization object, NOT their finite endpoint value.
    T must be an ordinary commutative Symbol. The exact coefficients use
    convergent alternating endpoints and no copied special-value tables.
    """
    if type(T) is not s.Symbol or T.is_commutative is not True:
        raise TypeError('T must be an ordinary commutative SymPy Symbol')
    word = _word(word)
    if len(word) > 12:
        raise UnsupportedExpressionError("Shuffle regularization exceeds the weight-12 resource budget")
    return _polynomial(word, T)


def shuffle_regularized_hpl_at_one(word):
    """Explicit constant term of the real endpoint polynomial; not a limit.

    In particular H_1 has regularized constant zero although it diverges.
    Existing convergent endpoint APIs retain their divergence errors. This
    helper is not used implicitly by the hypergeometric expansion engine.
    """
    T = s.Symbol('_regularization_T')
    return hpl_endpoint_polynomial(word, T).subs(T, 0)
