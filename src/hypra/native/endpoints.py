"""Unit endpoint identities and exact terminating hypergeometric sums.

Rewriting precedes Laurent planning because Gamma quotients and excess
recurrences can introduce regulator poles.
"""
from __future__ import annotations

from functools import cache

import sympy as s

from ..errors import UnsupportedExpressionError
from .planning import _positive


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
        return _positive(real_excess, assumptions)
    except (ValueError, TypeError):
        return False


def convergent_unit(hyper, epsilon):
    """Test strictly positive base excess for affine integer/half-integer pF(p-1)."""
    from .matched_half import half_parameters

    parsed = half_parameters(hyper, epsilon)
    if parsed is None or len(hyper.ap) != len(hyper.bq)+1:
        return False
    upper, lower = parsed
    if any(base.is_Integer and base <= 0 and slope == 0 for base, slope in lower):
        return False
    return sum(base for base, _ in lower) > sum(base for base, _ in upper)


def _raise_unit_excess(upper, lower, excess):
    """Sum the telescoping coefficient difference on the convergent side."""
    n = s.Dummy("index")
    denominator = s.prod(n+b for b in lower)
    remainder = s.Poly(s.expand(s.prod(n+a for a in upper) + (excess-n)*denominator), n)
    # The two leading degrees cancel. Express the remainder in falling
    # factorials so each term can be summed by reindexing a shifted pFq.
    terms = []
    for j in range(len(lower)-1, -1, -1):
        coefficient = remainder.nth(j)
        if coefficient == 0:
            continue
        remainder -= s.Poly(coefficient*s.ff(n, j), n)
        scale = s.prod(s.rf(a, j) for a in upper)/s.prod(s.rf(b, j+1) for b in lower)
        shifted = s.hyper([a+j for a in upper], [b+j+1 for b in lower], 1, evaluate=False)
        terms.append((coefficient*scale, shifted))
    return terms


def _reduce_unit_rank(upper, lower):
    """Remove an upper/lower pair whose ratio is polynomial in the index."""
    for a in upper:
        for b in lower:
            shift = a-b
            if not shift.is_Integer or not 1 <= shift <= 8:
                continue
            ap, bq = list(upper), list(lower)
            ap.remove(a)
            bq.remove(b)
            # Expand (b+n)_shift in the falling factorials of n.
            return s.Add(*(
                s.binomial(shift, j) * s.rf(b+j, shift-j)
                * s.prod(s.rf(value, j) for value in ap)
                / (s.rf(b, shift) * s.prod(s.rf(value, j) for value in bq))
                * s.hyper([value+j for value in ap], [value+j for value in bq],
                          1, evaluate=False)
                for j in range(int(shift)+1)
            ))
    return None


def _continue_unit(upper, lower, excess, epsilon):
    """Iterate the convergent identity as a meromorphic parameter identity."""
    base, slope = _affine(excess, epsilon)
    steps = int(s.floor(-base)) + 1
    if steps > 8 or len(upper) > 8:
        raise UnsupportedExpressionError("Unit continuation exceeds eight shifts or rank eight")
    if base.is_Integer and slope == 0:
        raise UnsupportedExpressionError("Unit continuation lies on a fixed excess pole")

    @cache
    def raise_until_convergent(atom):
        delta = s.expand(sum(atom.bq)-sum(atom.ap))
        if delta.subs(epsilon, 0) > 0:
            return {atom: s.S.One}
        result = {}
        for weight, shifted in _raise_unit_excess(atom.ap, atom.bq, delta):
            for leaf, coefficient in raise_until_convergent(shifted).items():
                result[leaf] = result.get(leaf, s.S.Zero) + weight*coefficient
        return {leaf: s.cancel(coefficient/delta) for leaf, coefficient in result.items()}

    atom = s.hyper(upper, lower, 1, evaluate=False)
    return s.Add(*(weight*leaf for leaf, weight in raise_until_convergent(atom).items()))


def rewrite_endpoints(expr, epsilon, assumptions=s.true, *, meromorphic=False, excess_shift=None):
    """Rewrite terminating sums and unit endpoints before Laurent planning.

    Other nonterminating arguments are unchanged. With meromorphic=True,
    continue supported unit endpoints in the parameters. A fixed nonpositive lower parameter is rejected before
    termination, since its definition can depend on a limiting convention.
    """
    if excess_shift is not None and (
        not meromorphic or type(excess_shift) is not int or excess_shift < 1
    ):
        raise ValueError("excess_shift requires meromorphic mode and a positive index")
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
        if excess_shift is not None and len(upper) == len(lower)+1 and lower:
            excess = s.expand(sum(lower)-sum(upper))
            if _fixed_nonpositive_integer(excess, epsilon):
                if excess_shift > len(lower):
                    raise UnsupportedExpressionError(
                        "Unit excess shift index exceeds the reduced lower parameter count"
                    )
                # Move transversely before any rank reduction or excess
                # raising. Lower parameters follow SymPy's canonical order.
                lower[excess_shift-1] += epsilon
                shifted = s.hyper(upper, lower, 1, evaluate=False)
                replacements[atom] = rewrite_endpoints(
                    shifted, epsilon, assumptions, meromorphic=True
                )
                continue
        if len(upper) == 1 and not lower:
            # Binomial theorem at its convergent endpoint; often obtained by
            # exact cancellation in a Gauss sum. This is not continuation.
            if meromorphic or _convergent(-upper[0], epsilon, assumptions):
                replacements[atom] = s.S.Zero
                continue
            raise UnsupportedExpressionError("Binomial z=1 convergence is not proven")
        if len(upper) != 2 or len(lower) != 1:
            if convergent_unit(reduced, epsilon):
                replacements[atom] = reduced
                continue
            if len(upper) == len(lower)+1:
                from .matched_half import half_parameters

                excess = s.expand(sum(lower)-sum(upper))
                parsed = half_parameters(reduced, epsilon)
                if parsed is not None and (meromorphic or _convergent(excess, epsilon, assumptions)):
                    lowered = _reduce_unit_rank(upper, lower) if meromorphic else None
                    replacements[atom] = (
                        rewrite_endpoints(lowered, epsilon, assumptions, meromorphic=True)
                        if lowered is not None else _continue_unit(upper, lower, excess, epsilon)
                    )
                    continue
            raise UnsupportedExpressionError(
                "Nonterminating z=1 requires supported parameters and proven positive excess"
            )
        a, b = upper
        c = lower[0]
        for parameter in (a, b, c):
            _affine(parameter, epsilon)
        excess = s.expand(c-a-b)
        if meromorphic and _fixed_nonpositive_integer(excess, epsilon):
            raise UnsupportedExpressionError("Gauss continuation lies on a fixed excess pole")
        if not meromorphic and not _convergent(excess, epsilon, assumptions):
            raise UnsupportedExpressionError("Gauss z=1 convergence is not proven on a regulator neighborhood or punctured side")
        replacements[atom] = s.gamma(c)*s.gamma(excess)/(s.gamma(c-a)*s.gamma(c-b))
    return expr.xreplace(replacements)
