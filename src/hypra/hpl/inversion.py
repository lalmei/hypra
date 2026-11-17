"""Explicit upper/lower inversion boundary values from real HPL transforms."""

from collections import defaultdict
from functools import lru_cache
from typing import cast

import sympy as s

from .._typing import Boundary, IntegerWord, Interval
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL
from .regularization import hpl_endpoint_polynomial
from .transforms import _validate

_ENDPOINT_LOG = s.Symbol("_inversion_endpoint_log")
_PULLBACK = {0: ((0, -1),), 1: ((0, 1), (1, 1)), -1: ((0, -1), (-1, 1))}


@lru_cache(maxsize=1024)
def _inversion_terms(word, sign):
    """Integrate the pulled-back suffix ODE and fix its endpoint constant.

    At z->1-, -log(1-(1/z +/- i0)) differs from -log(1-z)
    by +/- i*pi plus a vanishing log(z). The endpoint polynomial therefore
    determines the integration constant even for divergent leading-one words;
    no ordinary finite value is assigned to their singular endpoint.
    """
    if not word:
        return (((), s.S.One),)
    terms: defaultdict[tuple[int, ...], s.Expr] = defaultdict(lambda: s.S.Zero)
    for tail, coefficient in _inversion_terms(word[1:], sign):
        for letter, factor in _PULLBACK[word[0]]:
            terms[(letter, *tail)] += factor * coefficient
    endpoint = cast(
        s.Expr,
        hpl_endpoint_polynomial(word, _ENDPOINT_LOG).subs(_ENDPOINT_LOG, sign * s.I * s.pi),
    )
    constant = endpoint - sum(
        coefficient * hpl_endpoint_polynomial(w, _ENDPOINT_LOG).subs(_ENDPOINT_LOG, 0)
        for w, coefficient in terms.items()
    )
    terms[()] = s.expand(constant)
    return tuple((w, s.expand(c)) for w, c in terms.items() if c != 0)


def invert_hpl(
    word: IntegerWord,
    z: s.Symbol,
    *,
    boundary: Boundary,
    domain: Interval | None = None,
) -> s.Expr:
    """Return H_word(1/z +/- i0), for explicitly real 0<z<1.

    boundary='upper' selects +i0; 'lower' selects -i0. These choices refer
    to the argument 1/z, not the direction of the intermediate negative cut.
    Supports ternary words through weight eight, with exact endpoint constants,
    and repeated-letter words through weight 24.

    Args:
        word: Expanded ternary word; repeated-letter words admit higher weight.
        z: SymPy symbol for the interior coordinate.
        boundary: Side of approach to 1/z.
        domain: Declare `"0<z<1"` explicitly.
    """
    word = _validate(word, z, domain)
    if boundary not in ("upper", "lower"):
        raise UnsupportedExpressionError("Inversion boundary must be 'upper' or 'lower'")
    if not word:
        return s.S.One
    if len(set(word)) == 1:
        logarithm = HPL((0,), z)
        sign = 1 if boundary == "upper" else -1
        primitive = {
            0: -logarithm,
            1: logarithm + HPL((1,), z) + sign * s.I * s.pi,
            -1: HPL((-1,), z) - logarithm,
        }[word[0]]
        return primitive ** len(word) / s.factorial(len(word))
    if len(word) > 8:
        raise UnsupportedExpressionError("General inversion supports weight at most eight")
    sign = 1 if boundary == "upper" else -1
    return s.Add(
        *(
            coefficient * (HPL(term, z) if term else s.S.One)
            for term, coefficient in _inversion_terms(word, sign)
        )
    )
