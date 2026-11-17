"""Known-function HPL reductions derived from shuffle and path composition.

Nielsen functions retain their exact HPL convention. No upstream expression
table is copied; weight-four binary words reduce through their three Lyndon
generators, all of which have Nielsen form.
"""

from functools import lru_cache
from typing import cast

import sympy as s

from .._typing import Interval, as_integer_tuple
from ..errors import UnsupportedExpressionError
from ..nielsen import NielsenPolylog
from ..symbolic import HPL, reduce_known_hpl
from .classical import classical_hpl
from .notation import HPLPM, hpl_pm
from .reduction import lyndon_reduce


def _nielsen_word(word, argument):
    zeros = 0
    while zeros < len(word) and word[zeros] == 0:
        zeros += 1
    tail = word[zeros:]
    if zeros and tail and len(set(tail)) == 1 and tail[0] in (-1, 1):
        return tail[0] ** len(tail) * NielsenPolylog(zeros, len(tail), tail[0] * argument)
    return HPL(word, argument)


@lru_cache(maxsize=128)
def _binary_known(word, argument):
    """Shuffle-reduce a {0,1} word of weight at most four."""
    if not word:
        return s.S.One
    reduced = lyndon_reduce(HPL(word, argument))
    result = reduced.xreplace(
        {atom: _nielsen_word(tuple(atom.word), argument) for atom in reduced.atoms(HPL)}
    )
    result = reduce_known_hpl(cast(s.Expr, result))
    return result.xreplace(
        {
            atom: s.polylog(atom.n + 1, atom.argument)
            for atom in result.atoms(NielsenPolylog)
            if atom.p == 1
        }
    )


def _negative_binary_known(word, argument):
    # Relabel the shuffle alphabet only. The zero letter remains log(argument),
    # avoiding an unjustified log(-argument) substitution across its cut.
    positive_word = tuple(-letter for letter in word)
    reduced = lyndon_reduce(HPL(positive_word, argument))
    negative = reduced.xreplace(
        {
            atom: HPL(tuple(-letter for letter in as_integer_tuple(atom.word)), argument)
            for atom in reduced.atoms(HPL)
        }
    )
    result = negative.xreplace(
        {atom: _nielsen_word(tuple(atom.word), argument) for atom in negative.atoms(HPL)}
    )
    return reduce_known_hpl(cast(s.Expr, result))


def _opposite_known(word, argument):
    # u=(1+z)/2 sends the {-1,+1} kernels to the {0,1} kernels, but moves the
    # basepoint from zero to 1/2. Chen composition with the inverse first path
    # gives sum H_prefix(u)*(-1)^len(suffix)*H_reverse(suffix)(1/2).
    mapped = tuple(0 if letter == -1 else 1 for letter in word)
    return s.Add(
        *(
            (-1) ** (len(word) - split)
            * _binary_known(mapped[:split], (1 + argument) / 2)
            * _binary_known(mapped[split:][::-1], s.S.Half)
            for split in range(len(word) + 1)
        )
    )


def reduce_hpl_to_known_functions(expression: s.Expr, z: s.Symbol, *, domain: Interval) -> s.Expr:
    """Reduce the HPL 2 known-function families on the real interval 0<z<1.

    Covers all ternary words through weight three, every two-letter word
    through weight four, repeated letters, and signed Nielsen words 0^n(±1)^p.
    Higher mixed words remain HPLs. Outputs use logs, classical polylogarithms
    and NielsenPolylog with explicit principal continuations. Plus/minus HPLs
    are first expanded into their defining ordinary-letter linear combination.

    Args:
        expression: Exact expression containing ordinary or plus/minus HPLs.
        z: Common symbolic argument of the HPLs being reduced.
        domain: Explicit declaration `"0<z<1"`.
    """
    if not isinstance(expression, s.Expr) or not isinstance(z, s.Symbol):
        raise TypeError("Use a SymPy expression and symbolic argument")
    if domain != "0<z<1":
        raise UnsupportedExpressionError("Explicit domain='0<z<1' is required")
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("Known-function reduction requires exact finite input")
    ordinary = expression.xreplace(
        {atom: hpl_pm(atom.word, atom.argument) for atom in expression.atoms(HPLPM)}
    )
    result = classical_hpl(cast(s.Expr, ordinary), z, domain=domain)
    changes = {}
    for atom in result.atoms(HPL):
        if atom.argument != z:
            if atom.argument.has(z):
                raise UnsupportedExpressionError("First transform HPLs to the common argument")
            continue
        word = as_integer_tuple(atom.word)
        direct = _nielsen_word(word, z)
        if direct != atom:
            changes[atom] = direct
        elif len(word) <= 4:
            letters = set(word)
            if letters <= {0, 1}:
                changes[atom] = _binary_known(word, z)
            elif letters <= {0, -1}:
                changes[atom] = _negative_binary_known(word, z)
            elif letters <= {-1, 1}:
                changes[atom] = _opposite_known(word, z)
    return cast(s.Expr, result.xreplace(changes))
