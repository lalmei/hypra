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


def hpl_expression_endpoint_polynomial(expression, variable, T):
    """Map an exact polynomial in same-argument HPLs to its endpoint polynomial.

    ``variable -> 1-`` and ``T=-log(1-variable)``. Coefficients must be
    independent of both variables and are treated as fixed finite constants
    at the chosen specialization. No generic limit, integration interchange,
    argument transformation, or inverse-HPL asymptotics are performed.

    The input must be a SymPy Expr (not source text), polynomial in HPL atoms
    all at exactly ``variable``. Rational argument-dependent coefficients,
    even apparently removable ones, and unknown variable-dependent functions
    are outside this deliberately narrow ring. Products/sums are retained
    until polynomial substitution so divergent logarithms can cancel.
    """
    from .symbolic import HPL

    if (type(variable) is not s.Symbol or variable.is_commutative is not True
            or type(T) is not s.Symbol or T.is_commutative is not True or variable == T):
        raise TypeError('Use distinct ordinary commutative Symbols for variable and T')
    if not isinstance(expression, s.Expr):
        raise TypeError('expression must be an exact SymPy Expr, not source text')
    if expression.has(T, s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError('Input must be exact, finite, and independent of T')
    if s.count_ops(expression) > 4096:
        raise UnsupportedExpressionError('Endpoint expression exceeds the operation budget')
    atoms = sorted(expression.atoms(HPL), key=s.default_sort_key)
    if any(atom.argument != variable for atom in atoms):
        raise UnsupportedExpressionError('All HPLs must have exactly the endpoint variable as argument')
    if any(len(atom.word) > 12 for atom in atoms):
        raise UnsupportedExpressionError('Endpoint expression exceeds the weight-12 resource budget')
    generators = [s.Dummy('endpoint_hpl') for _ in atoms]
    replaced = expression.xreplace(dict(zip(atoms, generators, strict=True)))
    if replaced.is_commutative is not True:
        raise UnsupportedExpressionError('Endpoint polynomial requires commutative scalar coefficients')
    if replaced.has(variable):
        raise UnsupportedExpressionError('Coefficients must be independent of the endpoint variable')
    if not atoms:
        return expression
    # Bound expanded size before Poly expands products; this conservative bound
    # may reject a huge expression whose final terms happen to cancel.
    def term_budget(node):
        if not node.has(*generators):
            return 1
        if node in generators:
            return 1
        if node.is_Add:
            count = sum(term_budget(arg) for arg in node.args)
        elif node.is_Mul:
            count = 1
            for arg in node.args:
                count *= term_budget(arg)
                if count > 4096:
                    break
        elif node.is_Pow and node.exp.is_Integer and 0 <= node.exp <= 12:
            count = term_budget(node.base) ** int(node.exp)
        else:
            raise UnsupportedExpressionError('Input must be polynomial in HPLs, with no inverses')
        if count > 4096:
            raise UnsupportedExpressionError('Endpoint expression exceeds the expanded-term budget')
        return count
    term_budget(replaced)
    try:
        polynomial = s.Poly(replaced, *generators)
    except s.PolynomialError as exc:
        raise UnsupportedExpressionError('Input must be polynomial in HPLs') from exc
    for monomial, _ in polynomial.terms():
        if sum(degree * len(atom.word) for degree, atom in zip(monomial, atoms, strict=True)) > 12:
            raise UnsupportedExpressionError('Endpoint expression exceeds total HPL weight twelve')
    mapping = {generator: hpl_endpoint_polynomial(atom.word, T)
               for generator, atom in zip(generators, atoms, strict=True)}
    return s.expand(polynomial.as_expr().xreplace(mapping))
