"""Exact Nielsen generalized-polylogarithm notation, with explicit convention.

S_{n,p}(z) = H_{0^n,1^p}(z), n,p positive integers. This is a lossless
notation layer, not a general HPL-to-classical-polylogarithm reduction.
"""

import sympy as s
from sympy.core.function import ArgumentIndexError

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

    def __new__(cls, n, p, argument, **options):
        for index in (n, p):
            if isinstance(index, bool) or not isinstance(index, (int, s.Integer)) or index < 1:
                raise ValueError('Nielsen indices n and p must be exact positive integers')
        if n+p > 10000:
            raise ValueError('Nielsen weight exceeds the notation resource limit')
        return super().__new__(cls, s.Integer(n), s.Integer(p), argument, **options)

    @property
    def n(self):
        return self.args[0]

    @property
    def p(self):
        return self.args[1]

    @property
    def argument(self):
        return self.args[2]

    def as_hpl(self):
        """Losslessly expand into zeros followed by ones, without reduction."""
        return HPL((0,)*int(self.n)+(1,)*int(self.p), self.argument)

    def fdiff(self, argindex=3):
        if argindex != 3:
            raise ArgumentIndexError(self, argindex)
        if self.n > 1:
            return NielsenPolylog(self.n-1, self.p, self.argument)/self.argument
        return (-s.log(1-self.argument))**self.p/(s.factorial(self.p)*self.argument)

    def _eval_derivative(self, symbol):
        return self.fdiff(3)*s.diff(self.argument, symbol)

    def _eval_rewrite_as_HPL(self, n, p, argument, **hints):
        return self.as_hpl()

    def _eval_rewrite_as_polylog(self, n, p, argument, **hints):
        if p == 1:
            return s.polylog(n+1, argument)
        return self

    def _latex(self, printer):
        return (rf'S_{{{printer._print(self.n)},{printer._print(self.p)}}}'
                rf'\!\left({printer._print(self.argument)}\right)')


def _expression(expression):
    if type(expression) is int:
        expression = s.Integer(expression)
    if not isinstance(expression, s.Expr):
        raise TypeError('Provide an exact symbolic expression, not source text')
    return expression


def nielsen_to_hpl(expression):
    """Replace every Nielsen atom by its exact expanded HPL word."""
    expression = _expression(expression)
    return expression.xreplace({atom: atom.as_hpl() for atom in expression.atoms(NielsenPolylog)})


def hpl_to_nielsen(expression):
    """Replace only words 0^n 1^p (n,p positive); leave other HPLs untouched."""
    expression = _expression(expression)
    replacements = {}
    for atom in expression.atoms(HPL):
        word = tuple(atom.word)
        n = 0
        while n < len(word) and word[n] == 0:
            n += 1
        p = len(word)-n
        if n and p and all(letter == 1 for letter in word[n:]):
            replacements[atom] = NielsenPolylog(n, p, atom.argument)
    return expression.xreplace(replacements)


def reduce_nielsen(expression):
    """Opt-in p=1 reduction S_{n,1}(z)=Li_{n+1}(z); no other claims."""
    expression = _expression(expression)
    return expression.xreplace({atom: s.polylog(atom.n+1, atom.argument)
                               for atom in expression.atoms(NielsenPolylog) if atom.p == 1})
