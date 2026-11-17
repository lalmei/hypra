"""Exact origin-continuation HPL values at i, without fitted constant tables."""

from functools import lru_cache

import sympy as s

from .._typing import IntegerWord
from ..errors import UnsupportedExpressionError
from ..symbolic import HPL, reduce_known_hpl
from .classical import classical_hpl
from .multiple_polylog import MultiplePolylog
from .notation import HPLPM, PMWord, _pm_indices, letters_to_pm, pm_to_letters
from .numeric import _word
from .reduction import _extract


def _quarter_period_polylog(order, sign):
    """Split Li_n(±i) into its even eta and odd Dirichlet-beta sums."""
    if order == 2:
        beta = s.Catalan
    elif order % 2:
        half = (order - 1) // 2
        beta = (
            (-1) ** half
            * s.euler(2 * half)
            * s.pi**order
            / (4 ** (half + 1) * s.factorial(2 * half))
        )
    else:
        # Separate the two odd residue classes in the Hurwitz zeta. Its
        # relation to polygamma uses only the convergent defining series.
        beta = 2 * s.polygamma(order - 1, s.Rational(1, 4)) / (
            4**order * s.factorial(order - 1)
        ) - (1 - s.Rational(1, 2) ** order) * s.zeta(order)
    return -(1 - s.Rational(1, 2) ** (order - 1)) * s.zeta(order) / 2**order + sign * s.I * beta


def _reduce_at_i_constants(expression):
    expression = expression.xreplace(
        {
            atom: s.polylog(atom.args[0], s.simplify(atom.args[1]))
            for atom in expression.atoms(s.polylog)
        }
    )
    replacements = {}
    for atom in expression.atoms(s.polylog):
        order, argument = atom.args
        if order.is_Integer is not True or order < 2:
            continue
        order = int(order)
        if order == 3 and argument == (1 + s.I) / 2:
            # Compare the cubic antiderivative of H_011(i) with its exact
            # shuffle/Möbius/square reduction. Both integration constants are
            # fixed at the origin; their difference proves this conjugate
            # relation without fitting numerical values or assuming a basis.
            replacements[atom] = (
                -s.polylog(3, (1 - s.I) / 2)
                - 5 * s.pi**2 * s.log(2) / 96
                + s.log(2) ** 3 / 24
                + 35 * s.zeta(3) / 32
            )
        elif order == 3 and argument == s.S.Half:
            # The trilogarithm reflection identity at its midpoint.
            replacements[atom] = 7 * s.zeta(3) / 8 - s.pi**2 * s.log(2) / 12 + s.log(2) ** 3 / 6
        elif argument in (s.I, -s.I):
            replacements[atom] = _quarter_period_polylog(order, 1 if argument == s.I else -1)
        elif order == 2 and argument in ((1 + s.I) / 2, (1 - s.I) / 2):
            # Landen's identity Li2(z/(z-1))=-Li2(z)-log(1-z)^2/2,
            # at z=∓i. The radial path meets no transformed branch cut.
            original = -s.I if argument == (1 + s.I) / 2 else s.I
            replacements[atom] = (
                -_quarter_period_polylog(2, 1 if original == s.I else -1)
                - s.log(1 - original) ** 2 / 2
            )
        elif (s.Abs(argument) - 1).is_positive is True:
            # Jonquiere inversion on principal branches. Exact positive real
            # arguments take the same lower-cut boundary as SymPy polylog.
            logarithm = s.log(-argument)
            replacements[atom] = -((-1) ** order) * s.polylog(order, 1 / argument) - (
                (2 * s.pi * s.I) ** order
                / s.factorial(order)
                * s.bernoulli(order, s.S.Half + logarithm / (2 * s.pi * s.I))
            )
    result = expression.xreplace(replacements)
    # Inversion may introduce (1±i)/2 dilogarithms, so reduce those in a
    # second pass. Inversion strictly decreases the modulus and cannot cycle.
    if result != expression and result.has(s.polylog):
        result = _reduce_at_i_constants(result)
    return s.expand(result.xreplace({atom: s.expand_complex(atom) for atom in result.atoms(s.log)}))


@lru_cache(maxsize=512)
def _segment(word, start, displacement):
    """Iterated integral of a word with no singular final letter at start."""
    indices, arguments = [], []
    zeros, previous = 0, displacement
    for letter in word:
        pole = letter - start
        if pole == 0:
            zeros += 1
        else:
            indices.append(zeros + 1)
            arguments.append(s.expand(s.cancel(previous / pole)))
            zeros, previous = 0, pole
    sign = (-1) ** (len(indices) - word.count(1))
    return sign * MultiplePolylog(s.Tuple(*indices), s.Tuple(*arguments))


@lru_cache(maxsize=512)
def _first_half(word):
    if not word:
        return s.S.One
    midpoint = s.I / 2
    z = s.Symbol("_hpl_midpoint")
    extracted = _extract(word, z, False, True)
    return extracted.xreplace(
        {
            atom: s.log(midpoint) if atom.word == (0,) else _segment(tuple(atom.word), 0, midpoint)
            for atom in extracted.atoms(HPL)
        }
    )


def _at_i_multiple_polylogs(word):
    # Compose the first half with the reversal of the segment i -> i/2.
    # Both series have cumulative argument products with modulus <= 1/2.
    return _first_half(word) + s.Add(
        *(
            (-1) ** cut * _segment(word[:cut][::-1], s.I, -s.I / 2) * _first_half(word[cut:])
            for cut in range(1, len(word) + 1)
        )
    )


@lru_cache(maxsize=1093)
def _ordinary_at_i_raw(word: tuple[int, ...]) -> s.Expr:
    """Exact H_word(i) on continuation from the origin through the upper quadrant.

    All ternary words through weight six reduce to exact logarithms and
    classical or multiple polylogarithms. Pure repeated letters and simple
    classical words are also supported through weight 24.

    Args:
        word: Expanded ternary word evaluated at the positive imaginary unit.
    """
    if not word:
        return s.S.One
    z = s.Dummy("real_formula_argument", positive=True)
    if len(word) <= 3:
        # Obtain the proved real formula first, then continue it analytically.
        # Continue along z=i*t. All transformed arguments avoid their cuts,
        # except the cubic square-log identity: its argument 1-z**2 reaches
        # the lower polylog cut and log(z**2) the upper logarithm cut. Those
        # are exactly the principal boundary values used by SymPy here.
        # This does not widen classical_hpl's real-only public contract.
        formula = classical_hpl(HPL(word, z), z, domain="0<z<1")
    else:
        formula = reduce_known_hpl(HPL(word, z))
        if formula.has(HPL):
            if len(word) <= 6:
                return _at_i_multiple_polylogs(word)
            raise UnsupportedExpressionError(
                "At i, mixed words above weight six are not implemented"
            )
    value = formula.subs(z, s.I)
    return _reduce_at_i_constants(value)


@lru_cache(maxsize=1093)
def _ordinary_at_i_reduced(word: tuple[int, ...]) -> s.Expr:
    if len(word) > 6:
        return _ordinary_at_i_raw(word)
    from .at_i_reduction import at_i_weight_reductions

    basis = at_i_weight_reductions(len(word))
    return s.expand(
        s.Add(
            *(
                coefficient * basis[_pm_indices(pm_word)]
                for pm_word, coefficient in letters_to_pm(word).items()
            )
        )
    )


def hpl_at_i(word: IntegerWord) -> s.Expr:
    """Return the exact ordinary HPL at i in the HPL 2.0 constant basis.

    Through weight six, generated identities leave one, nine and 23 chosen
    HPLPM primitives at weights four, five and six. Their numerical evaluation
    uses independent convergent MPL expressions. Repeated single letters and
    simple classical words extend to weight 24.
    """
    return _ordinary_at_i_reduced(_word(word))


@lru_cache(maxsize=1093)
def _pm_at_i_raw(word: tuple[int, ...]) -> s.Expr:
    """Cached independent exact MPL expression, bypassing the canonical reducer.

    This private hook lets retained HPLPM at-i primitives evaluate numerically
    without calling their own symbolic reduction again.
    """
    terms = pm_to_letters(HPLPM.from_indices(word, s.I).word)
    return s.Add(
        *(coefficient * _ordinary_at_i_raw(letters) for letters, coefficient in terms.items())
    )


def hpl_pm_at_i(word: PMWord) -> s.Expr:
    """Return the exact plus/minus HPL at i in the canonical weight-six basis."""
    indices = _pm_indices(word)
    if len(word) > 6:
        raise UnsupportedExpressionError("At i, plus/minus words support weight at most six")
    from .at_i_reduction import at_i_weight_reductions

    return at_i_weight_reductions(len(indices))[indices]
