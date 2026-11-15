"""Exact real HPL transform t=(1-z)/(1+z), with endpoint-derived constants."""

from collections import defaultdict
from functools import lru_cache

import sympy as s

from .._typing import IntegerWord, Interval
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL
from .regularization import hpl_endpoint_polynomial
from .transforms import _validate

_MAX_MOBIUS_WEIGHT = 8
_ENDPOINT_LOG = s.Symbol("_mobius_endpoint_log")
_PULLBACK = {0: ((1, -1), (-1, -1)), 1: ((0, -1), (-1, 1)), -1: ((-1, -1),)}


@lru_cache(maxsize=512)
def _mobius_terms(word):
    if not word:
        return (((), s.S.One),)
    # T(M(z)) = -log(z) - log(2) + O(z). At the transformed origin
    # each nonempty HPL has zero constant term in its power/log expansion.
    constant = hpl_endpoint_polynomial(word, _ENDPOINT_LOG).subs(_ENDPOINT_LOG, -s.log(2))
    terms: defaultdict[tuple[int, ...], s.Expr] = defaultdict(lambda: s.S.Zero, {(): constant})
    for tail, coefficient in _mobius_terms(word[1:]):
        for letter, sign in _PULLBACK[word[0]]:
            terms[(letter, *tail)] += sign * coefficient
    return tuple((word, coefficient) for word, coefficient in terms.items() if coefficient != 0)


def mobius_hpl(word: IntegerWord, z: s.Symbol, *, domain: Interval | None = None) -> s.Expr:
    """Return H_word((1-z)/(1+z)) in HPLs at z on the real interval 0<z<1.

    Every expanded {-1,0,1} word through weight eight is supported.
    Endpoint constants are preserved
    as exact alternating multiple-zeta values when not reduced. The transform
    maps the real interval to itself; no complex continuation is claimed.

    Args:
        word: Expanded ternary word of weight at most eight.
        z: SymPy symbol for the interior coordinate.
        domain: Declare `"0<z<1"` explicitly.
    """
    word = _validate(word, z, domain)
    if len(word) > _MAX_MOBIUS_WEIGHT:
        raise UnsupportedExpressionError(
            f"Mobius transformation supports weight at most {_MAX_MOBIUS_WEIGHT}"
        )
    return s.Add(
        *(
            coefficient * (HPL(term, z) if term else s.S.One)
            for term, coefficient in _mobius_terms(word)
        )
    )
