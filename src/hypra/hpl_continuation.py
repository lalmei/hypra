"""Explicit real-axis boundary continuation, without arbitrary winding paths."""

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import _as_hpl
from .hpl_inversion import invert_hpl
from .hpl_transform import negate_hpl
from .symbolic import HPL


def analytic_continue_hpl(word, x, *, region, boundary):
    """Continue H_word(x +/- i0) in a caller-declared open real region.

    Regions are '0<x<1', '-1<x<0', 'x>1', 'x<-1'; boundary is mandatory
    'upper' (+i0) or 'lower' (-i0). This is a domain contract, not an inferred
    inequality. No value at x=0,+1,-1 or arbitrary winding is assigned.
    """
    if not isinstance(x, s.Symbol):
        raise TypeError("Continuation variable must be a SymPy Symbol")
    atom = HPL(word, x)
    word = tuple(map(int, atom.word))
    if region not in ("0<x<1", "-1<x<0", "x>1", "x<-1"):
        raise UnsupportedExpressionError("Unrecognized open real continuation region")
    if boundary not in ("upper", "lower"):
        raise UnsupportedExpressionError("Continuation boundary must be 'upper' or 'lower'")
    if len(word) > 24 or (region in ("x>1", "x<-1") and len(word) > 6):
        raise UnsupportedExpressionError("Continuation exceeds its supported weight limit")
    if not word:
        return s.S.One
    if region == "0<x<1":
        return atom
    t = s.Dummy("real_interval_argument", positive=True)
    if region == "x>1":
        return invert_hpl(word, t, boundary=boundary, domain="0<z<1").subs(t, 1 / x)
    negated = negate_hpl(word, t, domain="0<z<1")
    if boundary == "lower":
        negated = negated.xreplace({s.I: -s.I})
    if region == "-1<x<0":
        return negated.subs(t, -x)
    # x=-y reverses the boundary: upper x means lower y, and conversely.
    negated = _as_hpl(negated, t)
    for function in negated.atoms(s.Function):
        if function.has(t) and not isinstance(function, HPL):
            raise UnsupportedExpressionError("Unmapped variable function in continuation")
    atoms = negated.atoms(HPL)
    if any(
        term.xreplace({h: s.S.One for h in atoms}).has(t)
        for term in s.Add.make_args(s.expand(negated))
    ):
        raise UnsupportedExpressionError("Unmapped variable coefficient in continuation")
    opposite = "lower" if boundary == "upper" else "upper"
    return s.expand(
        negated.xreplace(
            {h: invert_hpl(h.word, t, boundary=opposite, domain="0<z<1") for h in atoms}
        )
    ).subs(t, -1 / x)
