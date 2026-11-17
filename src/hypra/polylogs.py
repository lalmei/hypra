"""Branch-explicit polylogarithm argument transformations from iterated integrals.

The six cross-ratio arguments in HypExp Appendix A.1 permute 0, 1 and infinity.
Their pulled-back one-forms therefore close on the binary HPL alphabet. Each
transformation below integrates those forms recursively and fixes its constant
using the endpoint logarithmic polynomial. No transformation table is used.
"""

from functools import lru_cache
from typing import cast

import sympy as s

from ._typing import Boundary, Interval
from .errors import UnsupportedExpressionError
from .hpl.regularization import hpl_endpoint_polynomial
from .nielsen import NielsenPolylog
from .symbolic import HPL, reduce_known_hpl

_T = s.Symbol("polylog_endpoint_log")
# Images of dx/x and dx/(1-x), expressed in dz/z and dz/(1-z).
_FORMS = {
    "reflection": {0: ((1, -1),), 1: ((0, -1),)},
    "fractional": {0: ((0, 1), (1, 1)), 1: ((1, -1),)},
    "inverse": {0: ((0, -1),), 1: ((0, 1), (1, 1))},
    "inverse_reflection": {0: ((1, 1),), 1: ((0, -1), (1, -1))},
    "reflection_inverse": {0: ((0, -1), (1, -1)), 1: ((0, 1),)},
}


@lru_cache(maxsize=2048)
def _word_transform(word: tuple[int, ...], mode: str, sign: int):
    """Return a linear HPL combination as (word, exact coefficient) pairs.

    Inversion and 1-1/z are normalized at z=1. The other maps are normalized
    at z=0. Endpoint finite parts fix integration constants without assigning
    finite values to divergent functions. For x>1 on the selected side,
    -log(1-x) has finite part sign*i*pi at the preimage of x=1.
    """
    if not word:
        return (((), s.S.One),)
    terms = {}
    for tail, coefficient in _word_transform(word[1:], mode, sign):
        for letter, factor in _FORMS[mode][word[0]]:
            key = (letter, *tail)
            terms[key] = terms.get(key, s.S.Zero) + factor * coefficient

    if mode in ("reflection", "inverse", "inverse_reflection"):
        phase = 0 if mode == "reflection" else sign * s.I * s.pi
        constant = hpl_endpoint_polynomial(word, _T).subs(_T, phase)
    else:
        # At x=0, only all-zero words have a logarithmic constant.
        constant = (
            (sign * s.I * s.pi) ** len(word) / s.factorial(len(word))
            if all(letter == 0 for letter in word)
            else s.S.Zero
        )
    if mode in ("inverse", "reflection_inverse"):
        constant -= s.Add(
            *(
                coefficient * hpl_endpoint_polynomial(key, _T).subs(_T, 0)
                for key, coefficient in terms.items()
            )
        )
    terms[()] = constant
    return tuple((key, s.expand(value)) for key, value in terms.items() if value != 0)


def transform_polylog_arguments(
    expression: s.Expr,
    z: s.Symbol,
    *,
    domain: Interval,
    boundary: Boundary | None = None,
) -> s.Expr:
    """Express classical and Nielsen polylogarithms at a common real argument.

    Supported arguments are z, 1-z, z/(z-1), 1-1/z, 1/z and 1/(1-z),
    with domain='0<z<1' and total weight at most eight per atom. Arguments
    above one require boundary='upper' or 'lower', referring to the argument
    of each original function. Negative arguments are analytic for these
    origin-regular words and need no side choice.

    The result uses exact binary HPLs at z, with simple words replaced by
    logarithms and classical polylogarithms. It covers the argument relations
    in HypExp Appendix A.1 without imposing its particular named-function
    basis or changing global evaluation rules. Constant arguments are retained;
    unknown variable arguments and unsupported orders raise explicitly.

    Args:
        expression: Exact expression containing SymPy polylog or NielsenPolylog.
        z: Common interior argument, a commutative SymPy Symbol.
        domain: Declare the real interval '0<z<1'.
        boundary: Optional upper or lower side for arguments above one.
    """
    if not isinstance(expression, s.Expr) or not isinstance(z, s.Symbol):
        raise TypeError("Use an exact SymPy expression and symbolic target argument")
    if z.is_commutative is not True:
        raise TypeError("The target argument must be commutative")
    if domain != "0<z<1":
        raise UnsupportedExpressionError("Explicit domain='0<z<1' is required")
    if boundary not in (None, "upper", "lower"):
        raise ValueError("boundary must be upper, lower, or omitted")
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("Polylogarithm transformation requires exact finite input")
    if s.count_ops(expression) > 4096:
        raise UnsupportedExpressionError("Polylogarithm expression exceeds the operation budget")
    candidates = (
        (z, "identity"),
        (1 - z, "reflection"),
        (z / (z - 1), "fractional"),
        (1 / z, "inverse"),
        (1 / (1 - z), "inverse_reflection"),
        (1 - 1 / z, "reflection_inverse"),
    )
    replacements = {}
    for atom in expression.atoms(s.polylog, NielsenPolylog):
        argument = atom.args[-1]
        if not argument.has(z):
            continue
        if isinstance(atom, NielsenPolylog):
            weight = int(atom.n + atom.p)
        else:
            order = atom.args[0]
            if not isinstance(order, s.Integer) or order < 1:
                raise UnsupportedExpressionError(
                    "Classical polylog orders must be positive integers"
                )
            weight = int(order)
        if weight > 8:
            raise UnsupportedExpressionError(
                "Polylogarithm transformation supports weight at most eight"
            )
        word = (
            (0,) * int(atom.n) + (1,) * int(atom.p)
            if isinstance(atom, NielsenPolylog)
            else (0,) * (weight - 1) + (1,)
        )
        mode = next((name for target, name in candidates if s.cancel(argument - target) == 0), None)
        if mode is None:
            raise UnsupportedExpressionError(f"Unsupported polylogarithm argument {argument}")
        if mode in ("inverse", "inverse_reflection") and boundary is None:
            raise UnsupportedExpressionError("Arguments above one require an explicit boundary")
        if mode == "identity":
            replacements[atom] = HPL(word, z)
        else:
            replacements[atom] = s.Add(
                *(
                    coefficient * (HPL(key, z) if key else s.S.One)
                    for key, coefficient in _word_transform(
                        word, mode, -1 if boundary == "lower" else 1
                    )
                )
            )
    return reduce_known_hpl(cast(s.Expr, expression.xreplace(replacements)))
