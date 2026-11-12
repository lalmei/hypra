"""Convergent endpoint identities and exact terminating hypergeometric sums.

No divergent sum is assigned its analytic continuation. Endpoint rewriting
precedes Laurent pole planning, since the resulting Gamma factors can have
regulator poles even when the original notation hides them.
"""
from __future__ import annotations

import sympy as s

from .errors import UnsupportedExpressionError


def _fixed_nonpositive_integer(value, epsilon):
    return not value.has(epsilon) and value.is_integer is True and value.is_nonpositive is True


def _affine(value, epsilon):
    base = value.subs(epsilon, 0)
    slope = s.diff(value, epsilon)
    if slope.has(epsilon) or s.expand(value-base-slope*epsilon) != 0:
        raise UnsupportedExpressionError("Gauss endpoint support requires affine regulator parameters")
    if base.has(s.oo, -s.oo, s.zoo, s.nan) or slope.has(s.oo, -s.oo, s.zoo, s.nan):
        raise UnsupportedExpressionError("Endpoint parameters must be finite")
    return base, slope


def _convergent(excess, epsilon, assumptions):
    base, slope = _affine(excess, epsilon)
    # A positive rational base gives an open neighborhood of convergence,
    # irrespective of the direction of a fixed finite complex slope.
    if isinstance(base, s.Rational) and base > 0:
        return True
    if base != 0 or slope == 0:
        return False
    # At the boundary, require an explicit punctured-side domain. Do not
    # infer convergence merely from the existence of a Gamma continuation.
    real_excess = s.re(excess)
    try:
        return s.ask(s.Q.positive(real_excess), assumptions) is True
    except (ValueError, TypeError):
        return False


def convergent_integer_unit(hyper, epsilon):
    """Strict positive integer-base excess for uniform unit-endpoint series.

    Moving lower poles are permitted meromorphically; fixed poles are not.
    This deliberately excludes zero base excess and all half-integer patterns.
    """
    from .integer import integer_parameters
    parsed = integer_parameters(hyper, epsilon)
    if parsed is None or len(hyper.ap) != len(hyper.bq)+1:
        return False
    upper, lower = parsed
    if any(base <= 0 and slope == 0 for base,slope in lower):
        return False
    return sum(base for base,_ in lower)-sum(base for base,_ in upper)>0


def rewrite_endpoints(expr, epsilon, assumptions=s.true):
    """Rewrite fixed terminating hypergeometric sums and convergent 2F1(1).

    Other nonterminating arguments are unchanged. Unsupported/divergent z=1
    atoms fail closed. A fixed nonpositive lower parameter is rejected before
    termination, since its definition can depend on a limiting convention.
    """
    if assumptions is None:
        assumptions = s.true
    if not isinstance(expr, s.Expr) or not isinstance(epsilon, s.Symbol):
        raise TypeError("Expected a SymPy expression and regulator Symbol")
    if not isinstance(assumptions, s.logic.boolalg.Boolean):
        raise TypeError("assumptions must be a SymPy Boolean")
    if s.simplify(assumptions) is s.false:
        raise UnsupportedExpressionError("Contradictory endpoint assumptions")
    replacements = {}
    for atom in expr.atoms(s.hyper):
        upper, lower, argument = list(atom.ap), list(atom.bq), atom.argument
        # Check before simplification or z=0 shortcuts: cancellation at a fixed
        # singular parameter selects a limiting convention we do not assume.
        if any(_fixed_nonpositive_integer(b, epsilon) for b in lower):
            raise UnsupportedExpressionError("Fixed nonpositive lower parameter has an ambiguous singular definition")
        for b in tuple(lower):
            if b in upper:
                upper.remove(b)
                lower.remove(b)
        reduced = s.hyper(upper, lower, argument, evaluate=False)
        terminating = [u for u in upper if _fixed_nonpositive_integer(u, epsilon)]
        if not terminating and argument != 1:
            if reduced != atom:
                replacements[atom] = reduced
            continue
        if atom.has(s.Float, s.oo, -s.oo, s.zoo, s.nan):
            raise UnsupportedExpressionError("Endpoint and terminating inputs must be exact and finite")
        if terminating:
            if any(not isinstance(u, s.Integer) for u in terminating):
                raise UnsupportedExpressionError("Termination degree must be an explicit integer")
            degree = min(-int(u) for u in terminating)
            if degree > 4096:
                raise UnsupportedExpressionError("Terminating sum exceeds the supported 4096-degree budget")
            replacements[atom] = s.Add(*(argument**n * s.prod(s.rf(a, n) for a in upper) /
                                         (s.factorial(n)*s.prod(s.rf(b, n) for b in lower))
                                         for n in range(degree+1)))
            continue
        if len(upper) == 1 and not lower:
            # Binomial theorem at its convergent endpoint; often obtained by
            # exact cancellation in a Gauss sum. This is not continuation.
            if _convergent(-upper[0], epsilon, assumptions):
                replacements[atom] = s.S.Zero
                continue
            raise UnsupportedExpressionError("Binomial z=1 convergence is not proven")
        if len(upper) != 2 or len(lower) != 1:
            if convergent_integer_unit(reduced, epsilon):
                replacements[atom] = reduced
                continue
            raise UnsupportedExpressionError("Nonterminating z=1 requires Gauss or integer parameters with strictly positive base excess")
        a, b = upper
        c = lower[0]
        for parameter in (a, b, c):
            _affine(parameter, epsilon)
        excess = s.expand(c-a-b)
        if not _convergent(excess, epsilon, assumptions):
            raise UnsupportedExpressionError("Gauss z=1 convergence is not proven on a regulator neighborhood or punctured side")
        replacements[atom] = s.gamma(c)*s.gamma(excess)/(s.gamma(c-a)*s.gamma(c-b))
    return expr.xreplace(replacements)
