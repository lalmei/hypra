"""Proved low-weight conversion to logarithms and classical polylogarithms."""
import sympy as s

from .errors import UnsupportedExpressionError
from .symbolic import HPL, reduce_known_hpl


def classical_hpl(expression, z, *, domain):
    """Reduce every weight-two HPL at z, plus standard all-weight simple words.

    Requires the declared real interval 0<z<1. Higher mixed words remain HPLs;
    this does not claim a complete high-weight classical-function basis.
    """
    if not isinstance(expression,s.Expr) or not isinstance(z,s.Symbol):
        raise TypeError("Use a SymPy expression and symbolic argument")
    if domain!='0<z<1':
        raise UnsupportedExpressionError("Explicit domain='0<z<1' is required")
    if expression.has(s.Float,s.nan,s.zoo,s.oo,-s.oo):
        raise UnsupportedExpressionError("Classical conversion requires exact finite input")
    result=reduce_known_hpl(expression)
    formulas={
        (1,0):-s.polylog(2,z)-s.log(z)*s.log(1-z),
        (-1,0):s.polylog(2,-z)+s.log(z)*s.log(1+z),
        (1,-1):-s.log(2)*s.log(1-z)+s.polylog(2,(1-z)/2)-s.polylog(2,s.S.Half),
        (-1,1):-s.log(2)*s.log(1+z)+s.polylog(2,(1+z)/2)-s.polylog(2,s.S.Half),
    }
    changes={}
    for atom in result.atoms(HPL):
        if atom.argument!=z:
            if atom.argument.has(z):
                raise UnsupportedExpressionError("First transform HPLs to the common argument")
            continue
        word=tuple(atom.word)
        if word in formulas:
            changes[atom]=formulas[word]
    return result.xreplace(changes)
