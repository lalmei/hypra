"""Select Python-native canonical hypergeometric kernels by parameter pattern."""
from __future__ import annotations

import sympy as s

from .colored_sums import canonical_colored_hyper_series
from .half_system import canonical_half_series
from .nested_sums import canonical_hyper_series


def canonical_bases(upper, lower):
    """Keep integer upper 1 where possible; unmatched half systems use 0."""
    r = sum(base.q == 2 for base, _ in upper)
    q = sum(base.q == 2 for base, _ in lower)
    if r == q:
        ap = [s.Rational(1, 2) if base.q == 2 else (s.S.Zero if base <= 0 else s.S.One)
              for base, _ in upper]
    else:
        ap = [s.Rational(1, 2) if base.q == 2 else s.S.Zero for base, _ in upper]
    return ap, [s.Rational(1, 2) if base.q == 2 else s.S.One for base, _ in lower]


def build_canonical(bases, slopes, z, eps, order):
    upper, lower = bases
    ua, lb = slopes
    r = sum(base == s.Rational(1, 2) for base in upper)
    q = sum(base == s.Rational(1, 2) for base in lower)
    if r == q == 0:
        return canonical_hyper_series(list(zip(upper, ua, strict=True)), lb, z, eps, order)
    if r == q:
        return canonical_colored_hyper_series(
            list(zip(upper, ua, strict=True)), list(zip(lower, lb, strict=True)), z, eps, order
        )
    if len(upper) == 2 and (r, q) == (2, 1):
        from .integer import _truncate
        a, b = ua
        c = lb[0]
        base = canonical_half_series([0, 0], [c-a, c-b], lower, lb, z, eps, order)
        factor = sum(((c-a-b)*eps*s.log(1-z))**k/s.factorial(k)
                     for k in range(order+1)) / s.sqrt(1-z)
        return _truncate(factor*base, eps, order)
    return canonical_half_series(upper, ua, lower, lb, z, eps, order)


def use_contiguous(upper, lower):
    """Select algebraic shifts where integral/tail reduction is expensive."""
    r = sum(base.q == 2 for base, _ in upper)
    q = sum(base.q == 2 for base, _ in lower)
    supported = r == q or (r, q) in ((0, 1), (1, 0)) or (
        len(upper) == 2 and (r, q) == (2, 1)
    )
    return supported and (r != q or any(base < 0 for base, _ in upper))
