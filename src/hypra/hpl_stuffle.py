"""Exact finite-cutoff quasi-shuffle algebra for convergent endpoint constants."""
from functools import lru_cache

import sympy as s

from .errors import UnsupportedExpressionError
from .symbolic import MZV, EulerSum

_MAX_WEIGHT = 12
_MAX_TERMS = 8192


@lru_cache(maxsize=4096)
def _stuffle(left, right):
    if not left:
        return ((right, 1),)
    if not right:
        return ((left, 1),)
    terms = {}
    merged = (left[0][0] + right[0][0], left[0][1] * right[0][1])
    for first, a, b in ((left[0], left[1:], right),
                        (right[0], left, right[1:]),
                        (merged, left[1:], right[1:])):
        for tail, count in _stuffle(a, b):
            word = (first,) + tail
            terms[word] = terms.get(word, 0) + count
    if len(terms) > _MAX_TERMS:
        raise UnsupportedExpressionError('Constant stuffle exceeds the term budget')
    return tuple(terms.items())


def stuffle_expand_constants(expression):
    """Linearize polynomials in MZV/EulerSum using exact strict stuffle products.

    Colors multiply when equal summation indices merge. The identity holds at
    finite common outer cutoff, then passes to the convergent outer-cutoff
    limits used by both constant classes, including alternating leading index
    one. No divergent or regularized constants are introduced. Scalar factors
    cannot themselves contain endpoint constants inside arbitrary functions.
    Output uses MZV for all-positive colors, EulerSum otherwise. Maximum total
    weight is twelve, powers eight, and intermediate linear terms 8192.
    """
    if not isinstance(expression, s.Expr):
        raise TypeError('Use an exact SymPy expression')
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError('Stuffle expansion requires exact finite input')

    def multiply(left, right):
        terms = {}
        for a, x in left.items():
            for b, y in right.items():
                if sum(i for i, _ in a + b) > _MAX_WEIGHT:
                    raise UnsupportedExpressionError('Constant stuffle exceeds weight twelve')
                for word, count in _stuffle(a, b):
                    terms[word] = terms.get(word, s.S.Zero) + x*y*count
                    if len(terms) > _MAX_TERMS:
                        raise UnsupportedExpressionError('Constant stuffle exceeds the term budget')
        return terms

    def polynomial(node):
        if isinstance(node, (MZV, EulerSum)):
            if sum(node.indices) > _MAX_WEIGHT:
                raise UnsupportedExpressionError('Constant stuffle exceeds weight twelve')
            colors = (1,)*len(node.indices) if isinstance(node, MZV) else node.colors
            return {tuple(zip(map(int, node.indices), map(int, colors))): s.S.One}
        if not node.has(MZV, EulerSum):
            return {(): node}
        if node.is_Add:
            result = {}
            for arg in node.args:
                for word, coefficient in polynomial(arg).items():
                    result[word] = result.get(word, s.S.Zero) + coefficient
                    if len(result) > _MAX_TERMS:
                        raise UnsupportedExpressionError('Constant stuffle exceeds the term budget')
            return result
        if node.is_Mul:
            result = {(): s.S.One}
            for arg in node.args:
                result = multiply(result, polynomial(arg))
            return result
        if node.is_Pow and node.exp.is_Integer and 0 <= node.exp <= 8:
            result, base = {(): s.S.One}, polynomial(node.base)
            for _ in range(int(node.exp)):
                result = multiply(result, base)
            return result
        raise UnsupportedExpressionError('Only polynomial endpoint-constant expressions are supported')

    def constant(word):
        if not word:
            return s.S.One
        indices, colors = zip(*word)
        return MZV(indices) if all(color == 1 for color in colors) else EulerSum(indices, colors)

    return s.Add(*(coefficient*constant(word) for word, coefficient in polynomial(expression).items()))
