"""Explicit upper/lower inversion boundary values from real HPL transforms."""

from typing import cast

import sympy as s

from .._typing import Boundary, IntegerWord, Interval
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL
from .algebra import _as_hpl
from .mobius import mobius_hpl
from .transforms import _validate, negate_hpl


def _map_hpl(expression, z, transform):
    expression = _as_hpl(expression, z)
    # Do not silently continue an unknown variable-dependent function.
    for atom in expression.atoms(s.Function):
        if atom.has(z) and not isinstance(atom, HPL):
            raise UnsupportedExpressionError(
                "Unrecognized variable function in inversion composition"
            )
    atoms = expression.atoms(HPL)
    if expression.xreplace({atom: s.Dummy() for atom in atoms}).has(z):
        raise UnsupportedExpressionError(
            "Unexpected variable-dependent coefficient in inversion composition"
        )
    return s.expand(
        expression.xreplace({atom: transform(atom.word) for atom in expression.atoms(HPL)})
    )


def invert_hpl(
    word: IntegerWord, z: s.Symbol, *, boundary: Boundary, domain: Interval | None = None,
) -> s.Expr:
    """Return H_word(1/z +/- i0), for explicitly real 0<z<1.

    boundary='upper' selects +i0; 'lower' selects -i0. These choices refer
    to the argument 1/z, not the direction of the intermediate negative cut.
    Supports ternary words through weight six, with exact endpoint constants.

    Args:
        word: Expanded ternary word of weight at most six.
        z: SymPy symbol for the interior coordinate.
        boundary: Side of approach to 1/z.
        domain: Declare `"0<z<1"` explicitly.
    """
    word = _validate(word, z, domain)
    if boundary not in ("upper", "lower"):
        raise UnsupportedExpressionError("Inversion boundary must be 'upper' or 'lower'")
    if len(word) > 6:
        raise UnsupportedExpressionError("Inversion supports weight at most six")
    expression = mobius_hpl(word, z, domain=domain)
    expression = _map_hpl(expression, z, lambda w: negate_hpl(w, z, domain=domain))
    lower = _map_hpl(expression, z, lambda w: mobius_hpl(w, z, domain=domain))
    # Every HPL(z), Euler sum, MZV and elementary real constant in this
    # representation is real on the declared interval. Only explicit I changes.
    return cast(s.Expr, lower if boundary == "lower" else lower.xreplace({s.I: -s.I}))
