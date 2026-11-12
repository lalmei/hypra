"""Proved low-weight conversion to logarithms and classical polylogarithms."""

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_reduce import extract_hpl_logs
from .symbolic import HPL, reduce_known_hpl


def classical_hpl(expression, z, *, domain):
    """Reduce all weight-two and binary weight-three HPLs at z.

    Requires the declared real interval 0<z<1. Higher mixed words remain HPLs;
    this does not claim a complete high-weight classical-function basis.
    """
    if not isinstance(expression, s.Expr) or not isinstance(z, s.Symbol):
        raise TypeError("Use a SymPy expression and symbolic argument")
    if domain != "0<z<1":
        raise UnsupportedExpressionError("Explicit domain='0<z<1' is required")
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("Classical conversion requires exact finite input")
    result = reduce_known_hpl(expression)
    formulas = {
        (1, 0): -s.polylog(2, z) - s.log(z) * s.log(1 - z),
        (-1, 0): s.polylog(2, -z) + s.log(z) * s.log(1 + z),
        (1, -1): -s.log(2) * s.log(1 - z) + s.polylog(2, (1 - z) / 2) - s.polylog(2, s.S.Half),
        (-1, 1): -s.log(2) * s.log(1 + z) + s.polylog(2, (1 + z) / 2) - s.polylog(2, s.S.Half),
    }
    # Binary weight three has only two irreducible words after extracting
    # leading ones and trailing zeros: 001=Li3(z), and 011 below. The latter
    # is fixed by derivative log(1-z)**2/(2*z) and its zero origin limit.
    binary_011 = (
        s.log(z) * s.log(1 - z) ** 2 / 2
        + s.log(1 - z) * s.polylog(2, 1 - z)
        - s.polylog(3, 1 - z)
        + s.zeta(3)
    )
    changes = {}
    for atom in result.atoms(HPL):
        if atom.argument != z:
            if atom.argument.has(z):
                raise UnsupportedExpressionError("First transform HPLs to the common argument")
            continue
        word = tuple(atom.word)
        if len(word) == 3 and all(letter in (0, 1) for letter in word):
            extracted = extract_hpl_logs(atom)
            changes[atom] = extracted.xreplace({HPL((0, 1, 1), z): binary_011})
        elif word in formulas:
            changes[atom] = formulas[word]
    return result.xreplace(changes)
