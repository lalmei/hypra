"""Multiple polylogarithm constants defined by absolutely convergent sums."""

from collections.abc import Iterable
from functools import lru_cache
from typing import cast

import mpmath as mp
import sympy as s
from sympy.core.evalf import prec_to_dps


@lru_cache(maxsize=512)
def _evaluate(indices, arguments, digits):
    ctx = mp.mp.clone()
    ctx.dps = digits + 16
    values = [
        ctx.mpc(*(str(part.evalf(ctx.dps)) for part in arg.as_real_imag())) for arg in arguments
    ]
    # Every cumulative argument product has modulus <= 1/2. For outer
    # index n the summands total at most n**(depth-1) * 2**(-n).
    depth = len(indices)
    cutoff = max(4 * depth, 4 * ctx.dps)
    while 4 * cutoff ** (depth - 1) * ctx.power(2, -cutoff) > ctx.power(10, -digits - 8):
        cutoff *= 2

    sums = [ctx.zero] * depth + [ctx.one]
    powers = [ctx.one] * depth
    for n in range(1, cutoff + 1):
        # Ascending j uses the previous cutoff for the inner strict sum.
        for j, (index, value) in enumerate(zip(indices, values, strict=True)):
            powers[j] *= value
            sums[j] += powers[j] * sums[j + 1] / n**index
    return ctx.nstr(sums[0].real, digits), ctx.nstr(sums[0].imag, digits)


class MultiplePolylog(s.Function):
    r"""Strict multiple polylogarithm Li_(indices)(arguments).

    The arguments here are exact Gaussian rationals with each cumulative
    product of modulus at most 1/2. Rewrite as Sum for its defining series.
    """

    nargs = 2

    @classmethod
    def eval(cls, indices, arguments):  # pyright: ignore[reportIncompatibleMethodOverride]
        if not isinstance(indices, s.Tuple) or not isinstance(arguments, s.Tuple):
            raise TypeError("Use SymPy tuples for multiple polylogarithm indices and arguments")
        if (
            not indices
            or len(indices) != len(arguments)
            or any(not isinstance(index, s.Integer) or index < 1 for index in indices)
        ):
            raise ValueError("Use one positive integer index per argument")
        cumulative = s.S.One
        for arg in arguments:
            if not isinstance(arg, s.Expr):
                raise TypeError("Multiple polylogarithm arguments must be Gaussian rationals")
            re, im = arg.as_real_imag()
            if not isinstance(re, s.Rational) or not isinstance(im, s.Rational):
                raise TypeError("Multiple polylogarithm arguments must be Gaussian rationals")
            cumulative = s.expand(cumulative * arg)
            re, im = cumulative.as_real_imag()
            if re**2 + im**2 > s.Rational(1, 4):
                raise ValueError("Cumulative argument products must have modulus at most 1/2")
        if len(indices) == 1:
            return s.polylog(indices[0], arguments[0])

    def _eval_rewrite_as_Sum(self, indices, arguments, **hints):
        variables = [s.Dummy("n", integer=True, positive=True) for _ in indices]
        terms = s.prod(
            arg**n / n**power for n, power, arg in zip(variables, indices, arguments, strict=True)
        )
        limits = [(variables[j], 1, variables[j - 1] - 1) for j in range(len(indices) - 1, 0, -1)]
        return s.Sum(terms, *limits, (variables[0], 1, s.oo))

    def _eval_evalf(self, prec):
        digits = prec_to_dps(prec) + 8
        indices = cast(Iterable[s.Integer], self.args[0])
        arguments = cast(Iterable[s.Expr], self.args[1])
        real, imag = _evaluate(tuple(map(int, indices)), tuple(arguments), digits)
        return s.Float(real, precision=prec) + s.I * s.Float(imag, precision=prec)
