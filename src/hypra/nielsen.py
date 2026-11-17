"""Exact Nielsen generalized-polylogarithm notation, with explicit convention.

S_{n,p}(z) = H_{0^n,1^p}(z), n,p positive integers. This is a lossless
notation layer, not a general HPL-to-classical-polylogarithm reduction.
"""

from typing import Self, cast

import sympy as s
from sympy.core.function import ArgumentIndexError

from ._typing import ExactInteger, ExactScalar
from .symbolic import HPL


class NielsenPolylog(s.Function):
    r"""S_{n,p}(z)=HPL((0,)*n+(1,)*p,z), n,p>=1.

    On 0<z<1 this equals
    (-1)**(n+p-1)/((n-1)! p!) * integral_0^1
    log(t)**(n-1) * log(1-z*t)**p dt/t.
    Elsewhere the symbol uses the matching principal HPL continuation, with
    no automatic branch transformations. Construction retains this notation.
    """

    nargs = 3

    def __new__(
        cls,
        n: ExactInteger,
        p: ExactInteger,
        argument: ExactScalar,
        **options: object,
    ) -> Self:
        """Construct the exact SymPy function.

        Args:
            n: Positive count of zero letters.
            p: Positive count of one letters.
            argument: SymPy argument of the generalized polylogarithm.
            **options: Keyword options forwarded to SymPy's function constructor.
        """
        for index in (n, p):
            if isinstance(index, bool) or not isinstance(index, (int, s.Integer)) or index < 1:
                raise ValueError("Nielsen indices n and p must be exact positive integers")
        if n + p > 10000:
            raise ValueError("Nielsen weight exceeds the notation resource limit")
        return cast(
            Self,
            super().__new__(cls, s.Integer(n), s.Integer(p), argument, **options),
        )

    @property
    def n(self) -> s.Integer:
        """Number of initial zeros in the expanded HPL word."""
        return cast(s.Integer, self.args[0])

    @property
    def p(self) -> s.Integer:
        """Number of final ones in the expanded HPL word."""
        return cast(s.Integer, self.args[1])

    @property
    def argument(self) -> s.Expr:
        """Argument of the generalized polylogarithm."""
        return cast(s.Expr, self.args[2])

    def as_hpl(self) -> HPL:
        """Losslessly expand into zeros followed by ones, without reduction."""
        return HPL((0,) * int(self.n) + (1,) * int(self.p), self.argument)

    def fdiff(self, argindex: int = 3) -> s.Expr:
        if argindex != 3:
            raise ArgumentIndexError(self, argindex)
        if self.n > 1:
            return NielsenPolylog(int(self.n) - 1, self.p, self.argument) / self.argument
        return (-s.log(1 - self.argument)) ** self.p / (s.factorial(self.p) * self.argument)

    def _eval_derivative(self, s):
        return self.fdiff(3) * self.argument.diff(s)

    def _eval_rewrite_as_HPL(self, n, p, argument, **hints):
        return self.as_hpl()

    def _eval_rewrite_as_polylog(self, n, p, argument, **hints):
        if p == 1:
            return s.polylog(n + 1, argument)
        return self

    def _eval_evalf(self, prec):
        """Evaluate using the principal HPL germ, retaining ambiguous cut values."""
        import mpmath as mp
        from sympy.core.evalf import prec_to_dps

        argument = self.argument
        if argument.is_number is not True or argument.is_finite is not True or self.n + self.p > 24:
            return self
        if argument == 0:
            return s.S.Zero
        if argument == 1:
            from .hpl.endpoints import hpl_at_one

            return hpl_at_one(self.as_hpl().word).evalf(prec_to_dps(prec) + 4)
        if argument.is_negative is True and argument != -1:
            from .hpl.boundary_numeric import evaluate_hpl_boundary

            digits = max(2, prec_to_dps(prec) + 8)
            if digits > 800:
                return self
            region = "-1<x<0" if (argument + 1).is_positive is True else "x<-1"
            result = evaluate_hpl_boundary(
                self.as_hpl().word, argument, region=region, boundary="upper", dps=digits
            )
            # Binary words ending in 1 are regular at zero; the negative ray
            # is not a cut for Nielsen functions. Their exact value is real.
            return s.Float(mp.nstr(result.value.real, n=digits), precision=prec)
        if argument == -1:
            reduced = reduce_nielsen(self)
            return reduced.evalf(prec_to_dps(prec) + 4) if reduced != self else self
        result = self.as_hpl()._eval_evalf(prec)
        return self if result.has(HPL) else result

    def _latex(self, printer):
        return (
            rf"S_{{{printer._print(self.n)},{printer._print(self.p)}}}"
            rf"\!\left({printer._print(self.argument)}\right)"
        )


def _expression(expression: ExactScalar) -> s.Expr:
    if type(expression) is int:
        expression = s.Integer(expression)
    if not isinstance(expression, s.Expr):
        raise TypeError("Provide an exact symbolic expression, not source text")
    return expression


def nielsen_to_hpl(expression: ExactScalar) -> s.Expr:
    """Replace every Nielsen atom by its exact expanded HPL word.

    Args:
        expression: Exact symbolic expression containing Nielsen atoms.
    """
    expression = _expression(expression)
    return cast(
        s.Expr,
        expression.xreplace({atom: atom.as_hpl() for atom in expression.atoms(NielsenPolylog)}),
    )


def hpl_to_nielsen(expression: ExactScalar) -> s.Expr:
    """Replace only words 0^n 1^p (n,p positive); leave other HPLs untouched.

    Args:
        expression: Exact symbolic expression containing HPL atoms.
    """
    expression = _expression(expression)
    replacements = {}
    for atom in expression.atoms(HPL):
        word = tuple(atom.word)
        n = 0
        while n < len(word) and word[n] == 0:
            n += 1
        p = len(word) - n
        if n and p and all(letter == 1 for letter in word[n:]):
            replacements[atom] = NielsenPolylog(n, p, atom.argument)
    return cast(s.Expr, expression.xreplace(replacements))


def reduce_nielsen(expression: ExactScalar) -> s.Expr:
    """Reduce p=1 to Li_{n+1}, and exact constants at zero and +/-1.

    Endpoint constants use proved MZV/Euler-sum reductions, retaining exact
    irreducible constants. Special-value reductions are bounded at weight 24.

    Args:
        expression: Exact symbolic expression containing Nielsen atoms.
    """
    expression = _expression(expression)
    from .hpl.alternating import hpl_alternating_at_one, reduce_euler_sum
    from .hpl.endpoints import hpl_at_one

    replacements = {}
    for atom in expression.atoms(NielsenPolylog):
        if atom.p == 1:
            replacements[atom] = s.polylog(atom.n + 1, atom.argument)
        elif atom.argument == 0:
            replacements[atom] = s.S.Zero
        elif atom.n + atom.p <= 24:
            if atom.argument == 1:
                replacements[atom] = hpl_at_one(atom.as_hpl().word)
            elif atom.argument == -1:
                word = (0,) * int(atom.n) + (-1,) * int(atom.p)
                replacements[atom] = (-1) ** atom.p * reduce_euler_sum(hpl_alternating_at_one(word))
    return cast(s.Expr, expression.xreplace(replacements))
