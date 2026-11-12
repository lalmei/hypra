"""Explicit expression-level dispatch for supported HPL argument maps."""
import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_inversion import invert_hpl
from .hpl_maps import fractional_hpl, square_hpl
from .hpl_mobius import mobius_hpl
from .hpl_transform import negate_hpl, reflect_hpl
from .symbolic import HPL


def transform_hpl_arguments(expression, z, *, domain, boundary=None):
    """Rewrite recognized HPL arguments in an expression to the common z.

    Requires domain='0<z<1'. Supported arguments: z, -z, 1-z, z²,
    (1-z)/(1+z), z/(z-1), 1/z. Cut-valued arguments require explicit
    boundary='upper' or 'lower'. Constant-argument HPLs remain unchanged.
    Unknown variable-dependent arguments fail rather than pretending to unify.
    """
    if not isinstance(expression,s.Expr) or not isinstance(z,s.Symbol):
        raise TypeError("Use an exact SymPy expression and symbolic target argument")
    if domain!='0<z<1':
        raise UnsupportedExpressionError("Explicit domain='0<z<1' is required")
    if boundary not in (None,'upper','lower'):
        raise ValueError("boundary must be upper, lower, or omitted")
    if expression.has(s.Float,s.nan,s.zoo,s.oo,-s.oo):
        raise UnsupportedExpressionError("Argument transformation requires exact finite input")
    changes={}
    for atom in expression.atoms(HPL):
        arg=atom.argument
        if arg==z or not arg.has(z):
            continue
        candidates=[(-z,'negate'),(1-z,'reflect'),(z**2,'square'),
                    ((1-z)/(1+z),'mobius'),(z/(z-1),'fractional'),(1/z,'inverse')]
        selected=next((name for target,name in candidates if s.cancel(arg-target)==0),None)
        if selected is None:
            raise UnsupportedExpressionError(f"Unsupported HPL argument {arg}")
        if selected in ('negate','fractional','inverse') and boundary is None:
            raise UnsupportedExpressionError("Cut-valued arguments require an explicit boundary")
        if selected=='negate':
            value=negate_hpl(atom.word,z,domain=domain)
            changes[atom]=value if boundary=='upper' else value.xreplace({s.I:-s.I})
        elif selected=='fractional':
            changes[atom]=fractional_hpl(atom.word,z,boundary=boundary,domain=domain)
        elif selected=='inverse':
            changes[atom]=invert_hpl(atom.word,z,boundary=boundary,domain=domain)
        else:
            function={'reflect':reflect_hpl,'square':square_hpl,'mobius':mobius_hpl}[selected]
            changes[atom]=function(atom.word,z,domain=domain)
    return expression.xreplace(changes)
