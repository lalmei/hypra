"""General 2F1 integer shifts of supported half-integer canonical systems."""

from __future__ import annotations

import sympy as s

from .errors import UnsupportedExpressionError
from .integer import _truncate
from .matched_half import half_parameters
from .planning import laurent_coefficient
from .symbolic import reduce_known_hpl


def _operations(upper, lower):
    bases = [s.Rational(1, 2) if base.q == 2 else s.S.Zero for base, _ in upper]
    c = s.Rational(1, 2) if lower[0][0].q == 2 else s.S.One
    operations = []
    for i, (target, _) in enumerate(upper):
        while bases[i] < target:
            loss = int(bases[i] == 0)
            operations.append(("upper_up", i, tuple(bases), c, loss))
            bases[i] += 1
        while bases[i] > target:
            loss = int(c - bases[i] == 0)
            operations.append(("upper_down", i, tuple(bases), c, loss))
            bases[i] -= 1
    target = lower[0][0]
    while c > target:
        operations.append(("lower_down", None, tuple(bases), c, int(c - 1 == 0)))
        c -= 1
    while c < target:
        loss = int(c - bases[0] == 0) + int(c - bases[1] == 0)
        operations.append(("lower_up", None, tuple(bases), c, loss))
        c += 1
    return operations


def extra_depth(h, eps):
    parsed = half_parameters(h, eps)
    if parsed is None or len(h.ap) != 2 or len(h.bq) != 1:
        return 0
    upper, lower = parsed
    return sum(op[-1] for op in _operations(upper, lower))


def expand_shifted_twof1(h, eps, order):
    """Canonical ODE systems plus algebraic contiguous parameter shifts.

    Independent formal slopes avoid illegal division by a zero slope. Apparent
    slope singularities are cancelled before the requested rational slopes are
    substituted. No parameter limit is numerically guessed.
    """
    from .half_system import canonical_half_series

    parsed = half_parameters(h, eps)
    if parsed is None or len(h.ap) != 2 or len(h.bq) != 1:
        return None
    upper, lower = parsed
    r = sum(base.q == 2 for base, _ in upper)
    q = sum(base.q == 2 for base, _ in lower)
    if (r, q) not in ((0, 1), (1, 0), (1, 1), (2, 1)):
        return None
    if lower[0][0].is_Integer and lower[0][0] <= 0 and lower[0][1] == 0:
        raise UnsupportedExpressionError("Fixed singular lower parameter")
    terminating = [
        -int(base) for base, slope in upper if base.is_Integer and base <= 0 and slope == 0
    ]
    if terminating:
        result = sum(
            h.argument**n
            / s.factorial(n)
            * s.prod(s.rf(a, n) for a in h.ap)
            / s.prod(s.rf(b, n) for b in h.bq)
            for n in range(min(terminating) + 1)
        )
        return _truncate(result, eps, order, int(lower[0][0].is_Integer and lower[0][0] <= 0))
    alpha, beta, gamma = s.symbols("hypra_alpha hypra_beta hypra_gamma")
    formal = (alpha, beta)
    z = s.Dummy("z", positive=True)
    canonical_upper = [s.Rational(1, 2) if base.q == 2 else s.S.Zero for base, _ in upper]
    canonical_lower = [s.Rational(1, 2) if lower[0][0].q == 2 else s.S.One]
    operations = _operations(upper, lower)
    depth = order + sum(op[-1] for op in operations)
    if (r, q) == (2, 1):
        base = canonical_half_series(
            [0, 0], [gamma - alpha, gamma - beta], [s.Rational(1, 2)], [gamma], z, eps, depth
        )
        prefactor = sum(
            ((gamma - alpha - beta) * eps * s.log(1 - z)) ** k / s.factorial(k)
            for k in range(depth + 1)
        ) / s.sqrt(1 - z)
        result = _truncate(prefactor * base, eps, depth)
    else:
        result = canonical_half_series(
            canonical_upper, formal, canonical_lower, [gamma], z, eps, depth
        )
    if result is None:
        return None
    theta = lambda f: z * s.diff(f, z)
    pole_bound = 0
    for operation, index, bases, cbase, loss in operations:
        a, b = [base + slope * eps for base, slope in zip(bases, formal, strict=True)]
        c = cbase + gamma * eps
        if operation == "upper_up":
            chosen = (a, b)[index]
            result = (theta(result) + chosen * result) / chosen
        elif operation == "upper_down":
            chosen, other = (a, b) if index == 0 else (b, a)
            result = ((c - chosen - other * z) * result + (1 - z) * theta(result)) / (c - chosen)
        elif operation == "lower_down":
            result = (theta(result) + (c - 1) * result) / (c - 1)
            if cbase - 1 == 0:
                pole_bound += 1
        else:
            result = c * ((1 - z) * theta(result) / z + (c - a - b) * result) / ((c - a) * (c - b))
        depth -= loss
        result = _truncate(result, eps, depth, pole_bound)
    actual = {alpha: upper[0][1], beta: upper[1][1], gamma: lower[0][1]}
    final = s.S.Zero
    for k in range(-pole_bound, order + 1):
        coefficient = s.cancel(laurent_coefficient(result, eps, k))
        value = coefficient.subs(actual, simultaneous=True)
        if value.has(s.nan, s.zoo, s.oo, -s.oo):
            raise UnsupportedExpressionError(
                "Degenerate contiguous slope limit was not resolved exactly"
            )
        final += value * eps**k
    return reduce_known_hpl(final).xreplace({z: h.argument})
