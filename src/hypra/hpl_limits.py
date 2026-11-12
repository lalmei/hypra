"""Combined finite z->1- limits using exact binary HPL reflection.

This is not termwise endpoint evaluation: rational/logarithmic singularities
are retained until all coefficients at the reflected origin have combined.
"""
import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import _as_hpl, _origin_constant
from .hpl_transform import reflect_hpl
from .symbolic import HPL


def hpl_limit_at_one(expression, z):
    """Prove the combined finite real endpoint, or raise on unsupported data.

    Coefficients must be rational in z, with ordinary HPLs/logs/polylogs at z.
    Reflection supports binary words of weight<=8. No complex branch choice or
    divergent regularization is inferred. A failure is not an infinite value.
    """
    if not isinstance(expression,s.Expr) or not isinstance(z,s.Symbol):
        raise TypeError('Expected a SymPy expression and symbolic argument')
    if expression.has(s.Float,s.nan,s.zoo,s.oo,-s.oo):
        raise UnsupportedExpressionError('Endpoint expression must be exact and finite')
    transformed=_as_hpl(expression,z)
    t=s.Dummy('endpoint_distance',positive=True)
    replacements={h:reflect_hpl(h.word,t,domain='0<z<1') for h in transformed.atoms(HPL)}
    reflected=transformed.xreplace(replacements).subs(z,1-t)
    return _origin_constant(reflected,t)
