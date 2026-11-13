"""Exact real near-one power/log expansions from endpoint constants and kernels."""

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_numeric import _word
from .hpl_regularized import hpl_endpoint_polynomial


def hpl_series_at_one(word, t, *, through=3):
    """Expand H_word(1-t) for real t->0+, including t**through.

    Uses T=-log(t) with exact polynomial coefficients. The explicit remainder
    is O(t**(through+1)*log(t)**weight), NOT an unjustified plain power Order.
    The input t must be an ordinary commutative Symbol; its positive approach
    is part of this API's domain contract. No complex branch is inferred.
    Weight <=8, through <=24, and a coefficient-operation budget are enforced.
    """
    if not isinstance(t, s.Symbol) or t.is_commutative is not True:
        raise TypeError('t must be an ordinary commutative SymPy Symbol')
    if type(through) is not int or not 0 <= through <= 24:
        raise UnsupportedExpressionError('Near-one order must be an integer from zero through 24')
    word = _word(word)
    if len(word) > 8:
        raise UnsupportedExpressionError('Near-one expansion supports weight at most eight')
    if not word:
        return s.S.One
    # A private ordinary Symbol is required by the endpoint API; its name is
    # distinct from t. All endpoint coefficients are parameter-free constants.
    T = s.Symbol('_near_one_log' if t.name != '_near_one_log' else '_near_one_log_')
    tail = [s.S.One] + [s.S.Zero]*through
    for offset in range(len(word)-1, -1, -1):
        head = word[offset]
        coefficients = [hpl_endpoint_polynomial(word[offset:], T)]
        for n in range(1, through+1):
            if head == 1:
                rhs = -tail[n]
            elif head == 0:
                rhs = -sum(tail[:n], s.S.Zero)
            else:
                rhs = -sum((tail[m]/s.Integer(2)**(n-m) for m in range(n)), s.S.Zero)
            # Invert n-D_T on the finite polynomial space, without integration
            # constants or divisions by symbolic parameters.
            value, derivative, power = s.S.Zero, s.expand(rhs), 1
            while derivative != 0:
                value += derivative/s.Integer(n)**power
                derivative = s.diff(derivative,T)
                power += 1
                if power > len(word)+2:
                    raise UnsupportedExpressionError('Near-one coefficient is not a finite log polynomial')
            value = s.expand(value)
            if s.count_ops(value) > 12000:
                raise UnsupportedExpressionError('Near-one coefficient operation budget exceeded')
            coefficients.append(value)
        tail = coefficients
    expression = sum(coefficient*t**n for n,coefficient in enumerate(tail))
    return s.expand(expression.subs(T,-s.log(t))) + s.Order(
        t**(through+1)*s.log(t)**len(word), (t,0))
