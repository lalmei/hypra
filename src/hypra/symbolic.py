"""Lossless harmonic-polylogarithm representation.

``HPL((a1, ..., an), z)`` uses the expanded Remiddi--Vermaseren alphabet
``ai in {-1, 0, 1}`` with kernels ``1/(1+t)``, ``1/t`` and ``1/(1-t)``.
In particular H_1(z) = -log(1-z), not log(1-z). The object deliberately
retains the symbolic function, including its word, without branch-changing
simplifications. Its analytic continuation is the principal continuation
from the real interval 0 < z < 1 unless an external calculation specifies
otherwise. It does not encode a separate continuation path.
"""

from __future__ import annotations

from typing import Self, cast

from sympy import (
    Basic,
    Dummy,
    Expr,
    Function,
    Integer,
    Order,
    Tuple,
    ceiling,
    factorial,
    log,
    polylog,
)
from sympy.core.function import ArgumentIndexError, PoleError

from ._typing import ExactScalar, IntegerWord


class HPL(Function):
    """A harmonic polylogarithm with an explicit, expanded letter word.

    Use :meth:`from_compressed` for HypExp/HPL's compressed integer indices:
    ``HPL.from_compressed((2, -1), z) == HPL((0, 1, -1), z)``.
    A compressed zero is one zero letter. No numerical or branch-dependent
    reduction is performed automatically.
    """

    nargs = 2

    def __new__(cls, word: IntegerWord, argument: ExactScalar, **options: object) -> Self:
        """Construct the exact SymPy function.

        Args:
            word: Expanded exact letters from {-1, 0, 1}, in integration order.
            argument: SymPy argument of the harmonic polylogarithm.
            **options: Keyword options forwarded to SymPy's function constructor.
        """
        # Validate before Function's cache: Python (1.0,) compares equal to
        # (1,), and otherwise a cached exact object can bypass eval validation.
        if not isinstance(word, (tuple, list, Tuple)):
            raise TypeError("HPL word must be a tuple/list of exact letters")
        if any(
            isinstance(letter, bool)
            or not isinstance(letter, (int, Integer))
            or letter not in (-1, 0, 1)
            for letter in word
        ):
            raise ValueError("HPL requires exact expanded letters -1, 0, or 1")
        return cast(Self, super().__new__(cls, Tuple(*word), argument, **options))

    # SymPy supplies a fixed set of positional arguments based on nargs.
    @classmethod
    def eval(cls, word, argument):  # pyright: ignore[reportIncompatibleMethodOverride]
        if not isinstance(word, Tuple):
            raise TypeError("HPL word must be a tuple of letters -1, 0, or 1")
        if any(not isinstance(letter, Integer) or letter not in (-1, 0, 1) for letter in word):
            raise ValueError("HPL requires expanded letters -1, 0, or 1; use from_compressed")

    @classmethod
    def from_compressed(cls, word: IntegerWord, argument: ExactScalar) -> Self:
        """Expand signed integer indices without changing their order.

            Args:
            word: Signed compressed indices. A standalone zero is one zero letter.
            argument: SymPy argument of the harmonic polylogarithm.

        Returns:
            An HPL with the corresponding expanded word.
        """
        if not isinstance(word, (tuple, list, Tuple)):
            raise TypeError("Compressed HPL indices must be a tuple or list of integers")
        letters = []
        for index in word:
            if isinstance(index, bool) or not isinstance(index, (int, Integer)):
                raise TypeError("Compressed HPL indices must be exact integers")
            index = int(index)
            if abs(index) > 10000:
                raise ValueError("Compressed HPL index exceeds the supported word length")
            if index == 0:
                letters.append(Integer(0))
            else:
                letters.extend([Integer(0)] * (abs(index) - 1))
                letters.append(Integer(1 if index > 0 else -1))
            if len(letters) > 10000:
                raise ValueError("HPL word exceeds the supported length")
        return cls(Tuple(*letters), argument)

    @property
    def word(self) -> Tuple:
        """Expanded letters, as a SymPy Tuple."""
        return cast(Tuple, self.args[0])

    @property
    def argument(self) -> Expr:
        """Argument of the harmonic polylogarithm."""
        return cast(Expr, self.args[1])

    def fdiff(self, argindex: int = 2) -> Expr:
        """Differentiate with respect to the argument, not the discrete word.

        The kernel is 1/z for a leading 0, 1/(1-z) for 1 and
        1/(1+z) for -1. The empty word is the constant function one.
        """
        if argindex != 2:
            raise ArgumentIndexError(self, argindex)
        if not self.word:
            return Integer(0)
        z = self.argument
        letter = cast(Integer, self.word[0])
        kernel = {Integer(0): 1 / z, Integer(1): 1 / (1 - z), Integer(-1): 1 / (1 + z)}[letter]
        tail = HPL(self.word[1:], z) if len(self.word) > 1 else Integer(1)
        return kernel * tail

    def _eval_derivative(self, s):
        # Function's generic implementation tries to differentiate Tuple(word),
        # but the letters are discrete constants. Only the argument varies.
        return self.fdiff(2) * self.argument.diff(s)

    def _eval_evalf(self, prec):
        """Evaluate interior real and off-cut complex principal values."""
        import mpmath as mp
        from sympy import Float, I
        from sympy.core.evalf import prec_to_dps

        z = self.argument
        if z.is_number is not True or z.is_finite is not True or len(self.word) > 24:
            return self
        digits = max(2, prec_to_dps(prec) + 8)
        word = self.word
        if z.is_real is False and digits <= 500:
            from .hpl.radial_continuation import continue_hpl

            value = continue_hpl(
                word, z, domain="principal-cut-plane", path="radial", dps=digits
            ).value
            return Float(mp.nstr(value.real, n=digits), precision=prec) + I * Float(
                mp.nstr(value.imag, n=digits), precision=prec
            )
        if z.is_positive is not True or (1 - z).is_positive is not True or digits > 1000:
            return self
        from .hpl.numeric import hpl_evalf

        value = hpl_evalf(word, z, dps=digits)
        return Float(mp.nstr(value, n=digits), precision=prec)

    def _zero_germ_argument(self, x, logx=None, cdir=0):
        """Return a leading argument and its positive rational valuation."""
        from .errors import UnsupportedExpressionError

        try:
            leading = self.argument.as_leading_term(x, logx=logx, cdir=cdir)
            coefficient, power = leading.as_coeff_exponent(x)
        except (AttributeError, NotImplementedError, ValueError, TypeError, PoleError) as exc:
            raise UnsupportedExpressionError(
                "Could not establish a zero-centered HPL argument germ"
            ) from exc
        if (
            coefficient.has(x)
            or coefficient.is_nonzero is not True
            or power.is_Rational is not True
            or power <= 0
        ):
            raise UnsupportedExpressionError(
                "HPL series requires an argument tending to zero with a known "
                "nonzero coefficient and positive rational power"
            )
        return leading, power

    def _eval_nseries(self, x, n, logx, cdir=0):
        """Exact zero-germ composition for analytic words and pure log words.

        SymPy's generic Function implementation differentiates the discrete
        Tuple word and cannot be used. The kernel-derived local series is used
        only where its omitted powers give a valid Order term after composition.
        """
        from .errors import UnsupportedExpressionError
        from .hpl.algebra import hpl_series

        if not self.word:
            return Integer(1)
        if not self.argument.has(x):
            return self
        if all(letter == 0 for letter in self.word):
            exact = log(self.argument) ** len(self.word) / factorial(len(self.word))
            return exact._eval_nseries(x, n, logx, cdir=cdir)
        center = self.argument.subs(x, 0)
        if (
            isinstance(center, Expr)
            and center.is_number is True
            and center.is_real is True
            and center.is_positive is True
            and (1 - center).is_positive is True
        ):
            from .hpl.taylor import hpl_taylor

            delta = self.argument - center
            leading = delta.as_leading_term(x, logx=logx, cdir=cdir)
            coefficient, power = leading.as_coeff_exponent(x)
            if (
                coefficient.has(x)
                or coefficient.is_nonzero is not True
                or not power.is_Rational
                or power <= 0
            ):
                raise UnsupportedExpressionError(
                    "Could not establish an interior Taylor argument germ"
                )
            degree = int(ceiling(max(n, 0) / power)) + 1
            local = Dummy("hpl_taylor")
            polynomial = hpl_taylor(self.word, center, local, through=degree).removeO()
            return polynomial.subs(local, self.argument)._eval_nseries(
                x, n, logx, cdir=cdir
            ) + Order(x**n, x)
        if self.word[-1] == 0:
            raise UnsupportedExpressionError(
                "Mixed trailing-zero HPL words have logarithmic germs; use "
                "hypra.hpl.algebra.hpl_series for explicit local terms"
            )
        _, power = self._zero_germ_argument(x, logx, cdir)
        # Expand strictly deeper than required: omitted analytic argument powers
        # have valuation > n, so appending O(x**n) cannot overstate accuracy.
        local_order = int(ceiling(max(n, 0) / power)) + 1
        local_variable = Dummy("hpl_germ", positive=True)
        germ = hpl_series(HPL(self.word, local_variable), local_variable, local_order)
        composed = germ.subs(local_variable, self.argument)
        return composed._eval_nseries(x, n, logx, cdir=cdir) + Order(x**n, x)

    # SymPy infers Self from its default hook, but accepts any leading expression.
    def _eval_as_leading_term(self, x, logx=None, cdir=0):  # pyright: ignore[reportIncompatibleMethodOverride]
        from .errors import UnsupportedExpressionError
        from .hpl.algebra import hpl_series

        if not self.word:
            return Integer(1)
        if not self.argument.has(x):
            return self
        if all(letter == 0 for letter in self.word):
            exact = log(self.argument) ** len(self.word) / factorial(len(self.word))
            return exact.as_leading_term(x, logx=logx, cdir=cdir)
        if self.word[-1] == 0:
            raise UnsupportedExpressionError(
                "Mixed trailing-zero HPL leading terms require explicit hpl_series"
            )
        leading, _ = self._zero_germ_argument(x, logx, cdir)
        # Each nonzero kernel increases the first analytic power by one; zero
        # kernels preserve it. Its coefficient comes from the exact germ engine.
        degree = sum(letter != 0 for letter in self.word)
        local_variable = Dummy("hpl_germ", positive=True)
        germ = hpl_series(HPL(self.word, local_variable), local_variable, degree + 1)
        return germ.as_leading_term(local_variable).subs(local_variable, leading)

    def _latex(self, printer):
        word = ",".join(printer._print(letter) for letter in self.word)
        return rf"H_{{{word}}}\!\left({printer._print(self.argument)}\right)"


def _constant_evalf(atom, prec, evaluator, arguments):
    """Explicit constant evalf, without importing numerical engines at startup."""
    import mpmath as mp
    from sympy import Float
    from sympy.core.evalf import prec_to_dps

    from .hpl.mobius import _MAX_MOBIUS_WEIGHT

    digits = max(2, prec_to_dps(prec) + 8)
    if sum(atom.indices) > _MAX_MOBIUS_WEIGHT or digits > 900:
        return atom
    value = evaluator(*arguments, dps=digits).value
    return Float(mp.nstr(value, n=digits), precision=prec)


class MZV(Function):
    """An exact convergent multiple zeta value in strict descending convention.

    MZV((s1,...,sd)) = sum_{n1>...>nd>=1} prod(nj**(-sj)), where every
    index is a positive integer and s1 >= 2. The empty word evaluates to one.
    No divergent regularization or numerical fitting is implied. Reduction to
    ordinary zeta values is opt-in through ``hypra.hpl.endpoints.reduce_mzv``.
    """

    nargs = 1

    def __new__(cls, indices: IntegerWord, **options: object) -> Expr:
        """Construct the exact SymPy function.

        Args:
            indices: Positive summation exponents with first index at least two.
            **options: Keyword options forwarded to SymPy's function constructor.
        """
        if not isinstance(indices, (tuple, list, Tuple)):
            raise TypeError("MZV indices must be a tuple/list of exact integers")
        if any(
            isinstance(index, bool) or not isinstance(index, (int, Integer)) or index <= 0
            for index in indices
        ):
            raise ValueError("MZV indices must be exact positive integers")
        if indices and indices[0] < 2:
            raise ValueError("A convergent MZV requires first index >= 2")
        return cast(Expr, super().__new__(cls, Tuple(*indices), **options))

    # SymPy supplies a fixed set of positional arguments based on nargs.
    @classmethod
    def eval(cls, indices):  # pyright: ignore[reportIncompatibleMethodOverride]
        if not isinstance(indices, Tuple):
            raise TypeError("MZV indices must be a tuple of positive integers")
        if not indices:
            return Integer(1)
        if any(not isinstance(index, Integer) or index <= 0 for index in indices):
            raise ValueError("MZV indices must be exact positive integers")
        if indices[0] < 2:
            raise ValueError("A convergent MZV requires first index >= 2")

    @property
    def indices(self) -> Tuple:
        """Strict descending summation exponents."""
        return cast(Tuple, self.args[0])

    def _eval_evalf(self, prec):
        from .hpl.constants import evaluate_mzv

        return _constant_evalf(self, prec, evaluate_mzv, (self.indices,))

    def _eval_derivative(self, s):
        return Integer(0)

    def _eval_is_real(self):
        return True

    def _eval_is_positive(self):
        return True

    def _latex(self, printer):
        indices = ",".join(printer._print(index) for index in self.indices)
        return rf"\zeta\!\left({indices}\right)"


class EulerSum(Function):
    r"""Strict colored sum sum_{n1>...>nd>=1} prod colors[j]**nj / nj**indices[j].

    Outer-cutoff convergence is used, including conditional convergence when
    the first index is one and its color is minus one. No regularization.
    Construction preserves indices and colors rather than reducing constants.
    """

    nargs = 2

    def __new__(cls, indices: IntegerWord, colors: IntegerWord, **options: object) -> Expr:
        """Construct the exact SymPy function.

        Args:
            indices: Positive summation exponents.
            colors: One sign (-1 or 1) per index.
            **options: Keyword options forwarded to SymPy's function constructor.
        """
        if not isinstance(indices, (tuple, list, Tuple)) or not isinstance(
            colors, (tuple, list, Tuple)
        ):
            raise TypeError("EulerSum indices and colors must be exact tuples/lists")
        if len(indices) != len(colors):
            raise ValueError("EulerSum requires one color per index")
        if any(isinstance(n, bool) or not isinstance(n, (int, Integer)) or n <= 0 for n in indices):
            raise ValueError("EulerSum indices must be positive exact integers")
        if any(
            isinstance(c, bool) or not isinstance(c, (int, Integer)) or c not in (-1, 1)
            for c in colors
        ):
            raise ValueError("EulerSum colors must be exact -1 or 1")
        if indices and indices[0] == 1 and colors[0] == 1:
            raise ValueError("Divergent EulerSum; no regularized value is assigned")
        return cast(Expr, super().__new__(cls, Tuple(*indices), Tuple(*colors), **options))

    # SymPy supplies a fixed set of positional arguments based on nargs.
    @classmethod
    def eval(cls, indices, colors):  # pyright: ignore[reportIncompatibleMethodOverride]
        if not indices:
            return Integer(1)

    @property
    def indices(self) -> Tuple:
        """Strict descending summation exponents."""
        return cast(Tuple, self.args[0])

    @property
    def colors(self) -> Tuple:
        """Sign raised to the corresponding summation index."""
        return cast(Tuple, self.args[1])

    def _eval_evalf(self, prec):
        from .hpl.constants import evaluate_euler_sum

        return _constant_evalf(self, prec, evaluate_euler_sum, (self.indices, self.colors))

    def _eval_derivative(self, s):
        return Integer(0)

    def _eval_is_real(self):
        return True

    def _latex(self, printer):
        indices = ",".join(printer._print(n) for n in self.indices)
        colors = ",".join(printer._print(c) for c in self.colors)
        return (
            rf"\operatorname{{ES}}\!\left[\begin{{matrix}}{indices}\\{colors}\end{{matrix}}\right]"
        )


def reduce_known_hpl(expression: ExactScalar) -> Expr:
    """Reduce selected HPL words to exact logarithms and classical polylogs.

    This is deliberately opt-in: construction/serialization retain HPL objects.
    Conventions are the principal branches continued from 0 < z < 1, with the
    same side of cuts in both HPL and output functions. This routine performs
    no argument transformations, no power/log simplification and no assertion
    that the expression is regular at a branch point.

    Supported identities are the empty word (one), all-zero, all-one,
    all-minus-one words, and zero strings followed by plus or minus one. Mixed words outside
    these patterns remain symbolic. Formulae hold away from singular endpoints
    and extend there only when the corresponding limit exists.
    """
    if type(expression) is int:
        expression = Integer(expression)
    if not isinstance(expression, Basic):
        raise TypeError("expression must be a SymPy expression")

    def reduce_one(atom):
        word, z = atom.word, atom.argument
        weight = len(word)
        if not word:
            return Integer(1)
        if all(letter == 0 for letter in word):
            return log(z) ** weight / factorial(weight)
        if all(letter == 1 for letter in word):
            return (-log(1 - z)) ** weight / factorial(weight)
        if all(letter == -1 for letter in word):
            return log(1 + z) ** weight / factorial(weight)
        if word[-1] in (-1, 1) and all(letter == 0 for letter in word[:-1]):
            return word[-1] * polylog(weight, word[-1] * z)
        return atom

    return cast(Expr, expression.replace(lambda node: isinstance(node, HPL), reduce_one))
