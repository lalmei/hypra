"""Exact origin-continuation HPL values at i, without fitted constant tables."""

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_classical import classical_hpl
from .hpl_numeric import _word
from .symbolic import HPL, reduce_known_hpl


def hpl_at_i(word):
    """Exact H_word(i) on continuation from the origin through the upper quadrant.

    All ternary words of weight <=2, pure repeated letters, and classical
    words 0**(n-1), +/-1 through weight24 are supported. Unsupported mixed
    words fail closed. Irreducible dilogarithms are retained exactly.
    """
    word = _word(word)
    if not word:
        return s.S.One
    z = s.Dummy("real_formula_argument", positive=True)
    if len(word) <= 2:
        # Obtain the proved real formula first, then continue it analytically.
        # Along z=i*t, 0<t<=1, every log/Li2 argument in these low-weight
        # formulas avoids its cut. This is not a complex-domain assertion
        # passed to classical_hpl's real-only public interface.
        formula = classical_hpl(HPL(word, z), z, domain="0<z<1")
    else:
        formula = reduce_known_hpl(HPL(word, z))
        if formula.has(HPL):
            raise UnsupportedExpressionError(
                "At i, mixed words above weight two are not implemented"
            )
    value = formula.subs(z, s.I)
    # The even/odd terms of Li2(±i) give -zeta(2)/8 ± i*Catalan.
    # Li2(1/2) follows from the reflection identity at 1/2.
    return s.expand(
        value.xreplace(
            {
                s.polylog(2, s.I): -(s.pi**2) / 48 + s.I * s.Catalan,
                s.polylog(2, -s.I): -(s.pi**2) / 48 - s.I * s.Catalan,
                s.polylog(2, s.S.Half): s.pi**2 / 12 - s.log(2) ** 2 / 2,
            }
        )
    )
