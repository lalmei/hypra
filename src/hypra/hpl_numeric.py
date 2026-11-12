"""High-precision explicitly branched HPL evaluation without nested numerical quadrature.

The default domain is real 0 < z < 1; explicit principal-disk mode permits
0 < |z| < 1 off the negative real axis, with principal logarithm. General words
use the convergent local logarithmic power series and an analytic tail bound.
Repeated-letter and classical-polylogarithm words use established mpmath
functions. Every result is recomputed at a higher working precision to check
stability; neither that check nor the tail bound is an interval-arithmetic
certificate of all roundoff errors.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import mpmath as mp
import sympy as s

from .errors import UnsupportedExpressionError

_DECIMAL = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z")
_MAX_WEIGHT = 24


class HPLNumericalConvergenceError(ArithmeticError):
    """The requested precision was not reached within the numerical budget."""


@dataclass(frozen=True)
class HPLNumericResult:
    """Numerical value and honest truncation/precision metadata.

    ``tail_bound`` bounds mathematical series truncation (zero for direct
    elementary/classical calls), not floating-point roundoff. ``stability_gap``
    is the difference between two working-precision calculations, not a proof
    of correctness. Requested precision cannot recover digits absent in input.
    """

    value: object
    tail_bound: object
    stability_gap: object
    terms: int
    dps: int
    working_dps: int
    method: str
    domain: str = "real 0 < z < 1; principal real HPL branch"


def _word(value):
    if not isinstance(value, (tuple, list, s.Tuple)):
        raise TypeError("Numerical HPL input must be a word of exact -1, 0, 1 letters")
    if any(
        isinstance(letter, bool)
        or not isinstance(letter, (int, s.Integer))
        or letter not in (-1, 0, 1)
        for letter in value
    ):
        raise UnsupportedExpressionError("Numerical HPL words require exact -1, 0, 1 letters")
    if len(value) > _MAX_WEIGHT:
        raise UnsupportedExpressionError(f"Numerical HPL weight exceeds {_MAX_WEIGHT}")
    return tuple(int(letter) for letter in value)


def _argument(value, ctx, domain):
    if isinstance(value, bool) or (domain == "real" and isinstance(value, (complex, mp.mpc))):
        raise TypeError("Numerical HPL argument must be real")
    if isinstance(value, s.Rational):
        result = ctx.mpf(int(value.p)) / int(value.q)
    elif isinstance(value, s.Expr):
        if (value.free_symbols or value.is_finite is not True
                or (domain == "real" and value.is_real is not True)):
            raise UnsupportedExpressionError(
                "Numerical HPL argument must be finite and satisfy the explicit domain"
            )
        real, imag = value.as_real_imag()
        result = ctx.mpc(
            ctx.mpf(str(s.N(real, ctx.dps + 10))),
            ctx.mpf(str(s.N(imag, ctx.dps + 10))),
        ) if domain == "principal-disk" else ctx.mpf(str(s.N(real, ctx.dps + 10)))
    elif isinstance(value, str):
        if len(value) > 4096 or not _DECIMAL.fullmatch(value):
            raise ValueError("Use a decimal string, not executable text, for numerical input")
        result = ctx.mpf(value)
    elif domain == "principal-disk" and isinstance(value, (complex, mp.mpc)):
        result = ctx.mpc(value)
    elif isinstance(value, (int, float, mp.mpf)):
        result = ctx.mpf(value)
    else:
        raise TypeError("Use an exact SymPy number, decimal string, integer, float or mpmath real")
    valid = (0 < result < 1) if domain == "real" else (
        0 < abs(result) < 1 and not (ctx.im(result) == 0 and ctx.re(result) <= 0)
    )
    if not ctx.isfinite(result) or not valid:
        raise UnsupportedExpressionError(
            "Numerical HPL requires real 0 < z < 1, or explicit principal-disk "
            "mode with 0 < |z| < 1 off the negative real axis"
        )
    return result


def _integrate_polynomial(coefficients, power, ctx):
    """Q with (power + d/dL) Q = P, represented by ascending powers of L."""
    result = [ctx.zero] * len(coefficients)
    following = ctx.zero
    for degree in range(len(coefficients) - 1, -1, -1):
        following = (coefficients[degree] - (degree + 1) * following) / power
        result[degree] = following
    return result


def _polynomial_value(coefficients, logarithm, ctx):
    value = ctx.zero
    for coefficient in reversed(coefficients):
        value = value * logarithm + coefficient
    return value


def _series(word, z, target_digits, max_terms, ctx):
    """O(weight * log-degree * terms) recursion, with constant-in-terms memory."""
    weight = len(word)
    trailing = 0
    for letter in reversed(word):
        if letter != 0:
            break
        trailing += 1
    logarithm = ctx.log(z)
    # All coefficient polynomials have degree at most the trailing-zero count.
    # Starting from L^d/d!, each integration has coefficient norm <= C_d/n,
    # while a nonzero kernel sums n preceding coefficients. Consequently the
    # uniform norm bound B below gives tail <= B*|z|^(N+1)/(1-|z|).
    c_degree = sum(
        math.factorial(trailing) // math.factorial(trailing - j) for j in range(trailing + 1)
    )
    coefficient_bound = (
        ctx.mpf(c_degree) ** (weight - trailing)
        * max(ctx.one, abs(logarithm)) ** trailing
        / math.factorial(trailing)
    )
    tolerance = ctx.power(10, -target_digits)
    size = trailing + 1
    # Each suffix has its previous coefficient polynomial. Pure-zero suffixes
    # contribute only n=0 and are initialized exactly, including empty=1.
    previous = []
    for start in range(weight + 1):
        coefficients = [ctx.zero] * size
        suffix = word[start:]
        if all(letter == 0 for letter in suffix):
            coefficients[len(suffix)] = ctx.one / math.factorial(len(suffix))
        previous.append(coefficients)
    cumulative = [[ctx.zero] * size for _ in range(weight)]
    total, z_power = ctx.zero, ctx.one
    for n in range(1, max_terms + 1):
        current = [[ctx.zero] * size for _ in range(weight + 1)]
        for start in range(weight - 1, -1, -1):
            if start >= weight - trailing:
                continue  # Pure-zero suffixes have no positive-power coefficients.
            head = word[start]
            if head == 0:
                integrand = current[start + 1]
            else:
                cumulative[start] = [
                    previous[start + 1][degree] + head * cumulative[start][degree]
                    for degree in range(size)
                ]
                integrand = cumulative[start]
            current[start] = _integrate_polynomial(integrand, n, ctx)
        z_power *= z
        total += z_power * _polynomial_value(current[0], logarithm, ctx)
        previous = current
        # A factor two leaves numerical slack when evaluating the exact analytic
        # majorant at finite working precision; this is not an interval bound.
        tail = 2 * coefficient_bound * abs(z_power) * abs(z) / (1 - abs(z))
        if total != 0 and tail <= tolerance * abs(total):
            return total, tail, n, "logarithmic-power-series"
    raise HPLNumericalConvergenceError(
        f"HPL series did not meet the relative truncation target after {max_terms} terms; "
        "increase max_terms or use an argument farther from the unit circle"
    )


def _evaluate(word, argument, dps, target_digits, max_terms, domain):
    # Clone avoids modifying the caller's global mpmath precision.
    ctx = mp.mp.clone()
    ctx.dps = dps
    z = _argument(argument, ctx, domain)
    weight = len(word)
    if not word:
        return ctx.one, ctx.zero, 0, "empty-word"
    if all(letter == 0 for letter in word):
        return ctx.log(z) ** weight / math.factorial(weight), ctx.zero, 0, "logarithm"
    if all(letter == 1 for letter in word):
        return (-ctx.log1p(-z)) ** weight / math.factorial(weight), ctx.zero, 0, "logarithm"
    if all(letter == -1 for letter in word):
        return ctx.log1p(z) ** weight / math.factorial(weight), ctx.zero, 0, "logarithm"
    if all(letter == 0 for letter in word[:-1]):
        sign = word[-1]
        return sign * ctx.polylog(weight, sign * z), ctx.zero, 0, "classical-polylogarithm"
    return _series(word, z, target_digits, max_terms, ctx)


def evaluate_hpl(word, argument, *, dps=50, max_terms=100000, domain="real"):
    """Evaluate a finite {0,±1} word in an explicitly selected branch domain.

    The requested ``dps`` is checked by a stricter analytic truncation target
    and agreement of calculations at dps+20 and dps+40. This is high-precision
    numerical evidence, not rigorous interval certification. Finite-precision
    float/mpmath inputs are evaluated at their stored value; prefer an exact
    SymPy Rational or decimal string to specify the intended point precisely.

    The default domain="real" requires 0 < z < 1. Explicit
    domain="principal-disk" allows 0 < |z| < 1 off the negative real axis,
    using principal log(z). No continuation beyond this slit disk is applied.
    Endpoints are rejected; use hpl_at_one for supported exact endpoints.
    """
    word = _word(word)
    if domain not in ("real", "principal-disk"):
        raise ValueError("domain must be real or principal-disk")
    if type(dps) is not int or not 2 <= dps <= 1000:
        raise ValueError("dps must be an integer between 2 and 1000")
    if type(max_terms) is not int or not 1 <= max_terms <= 1000000:
        raise ValueError("max_terms must be an integer between 1 and 1000000")
    first, _, _, _ = _evaluate(word, argument, dps + 20, dps + 8, max_terms, domain)
    second, tail, terms, method = _evaluate(word, argument, dps + 40, dps + 8, max_terms, domain)
    ctx = mp.mp.clone()
    ctx.dps = dps + 40
    gap = abs(ctx.convert(first) - second)
    if gap > ctx.power(10, -dps - 4) * abs(second):
        raise HPLNumericalConvergenceError(
            "HPL numerical result was not stable under increased working precision"
        )
    domain_description = (
        "real 0 < z < 1; principal real HPL branch" if domain == "real" else
        "principal slit disk 0 < |z| < 1, -pi < arg(z) < pi; principal log"
    )
    return HPLNumericResult(second, tail, gap, terms, dps, dps + 40, method, domain_description)


def hpl_evalf(word, argument, *, dps=50, max_terms=100000, domain="real"):
    """Return only the mpmath value; see evaluate_hpl for diagnostics."""
    return evaluate_hpl(word, argument, dps=dps, max_terms=max_terms, domain=domain).value
