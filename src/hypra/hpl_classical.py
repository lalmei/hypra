"""Proved low-weight conversion to logarithms and classical polylogarithms."""

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_reduce import _extract, extract_hpl_logs
from .symbolic import HPL, reduce_known_hpl


def classical_hpl(expression, z, *, domain):
    """Reduce weight two and all weight-three words using at most two letters.

    Requires the declared real interval 0<z<1. Weight-three words using all three letters and higher mixed words remain HPLs;
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
    # A real formula from u=1/(1+z), rather than substitution across the
    # Li_n(1+z) cuts in the positive-letter formula.
    u = 1 / (1 + z)
    L = s.log(1 + z)
    negative_011 = (
        s.zeta(3) - s.polylog(3, u) - L * s.polylog(2, u) + L**2 * s.log(z) / 2 - L**3 / 3
    )

    # The real substitution u=(1∓z)/2 evaluates H_{∓1,±1,±1}.
    def opposite_repeat(sign):
        logarithm = s.log(1 - sign * z)
        argument = (1 - sign * z) / 2
        return (
            s.polylog(3, s.S.Half)
            + logarithm**2 * s.log((1 + sign * z) / 2) / 2
            + logarithm * s.polylog(2, argument)
            - s.polylog(3, argument)
        )

    minus_plus_plus = opposite_repeat(1)
    plus_minus_minus = -opposite_repeat(-1)
    minus_minus_plus = (
        -(s.log(1 + z) ** 2) * s.log(1 - z) / 2
        - s.log(1 + z) * formulas[(1, -1)]
        + plus_minus_minus
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
        elif len(word) == 3 and all(letter in (0, -1) for letter in word):
            # Relabel only the abstract shuffle alphabet, not the argument.
            # This extracts leading -1 and trailing zero factors without any
            # complex negation or logarithmic branch substitution.
            binary = tuple(-letter for letter in word)
            extracted = _extract(binary, z, True, True)
            relabeled = extracted.xreplace(
                {h: HPL(tuple(-letter for letter in h.word), z) for h in extracted.atoms(HPL)}
            )
            changes[atom] = reduce_known_hpl(relabeled).xreplace(
                {HPL((0, -1, -1), z): negative_011}
            )
        elif len(word) == 3 and all(letter in (-1, 1) for letter in word):
            # Abstract two-letter shuffle extraction, not argument relabeling.
            binary = tuple(0 if letter == -1 else 1 for letter in word)
            extracted = _extract(binary, z, True, True)
            relabeled = extracted.xreplace(
                {
                    h: HPL(tuple(-1 if letter == 0 else 1 for letter in h.word), z)
                    for h in extracted.atoms(HPL)
                }
            )
            replacements = {HPL(w, z): f for w, f in formulas.items()}
            replacements.update(
                {HPL((-1, -1, 1), z): minus_minus_plus, HPL((-1, 1, 1), z): minus_plus_plus}
            )
            changes[atom] = reduce_known_hpl(relabeled.xreplace(replacements))
        elif word in formulas:
            changes[atom] = formulas[word]
    return result.xreplace(changes)
