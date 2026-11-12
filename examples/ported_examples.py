"""Independent exact benchmarks for native HypExp-style parameter classes.

These formulas are derived from the defining Pochhammer series (DLMF16.2.1),
not copied from HypExp code, tables or notebook output. HypExp (hep-ph/0507094)
and HypExp2 (0708.2443) motivate integer/half-integer expansion classes; no
claim is made that the papers print these precise examples. Official snapshot
installer tests only test package path variables, not mathematical results.

The geometric argument factor gives uniform convergence on |z|<=r<1 and
small regulator disks, justifying the finite-order coefficient expansion
under these sums. The endpoint example instead uses the convergent Gauss
identity on epsilon<0; no endpoint/integral interchange is assumed.

Sources: https://dlmf.nist.gov/16.2.E1 ; https://dlmf.nist.gov/15.4.E20 ;
https://arxiv.org/abs/hep-ph/0507094 ; https://arxiv.org/abs/0708.2443 .

Run: python examples/ported_examples.py
"""
from dataclasses import dataclass

import sympy as s

from hypra import HPL, NativeConfig, expand_epsilon
from hypra.hpl_reduce import extract_hpl_logs

EPSILON = s.Symbol('epsilon', real=True)
Z = s.Symbol('z', positive=True)
DOMAIN = s.And(Z>0,Z<1)


@dataclass(frozen=True)
class Example:
    name: str
    expression: s.Expr
    expected: s.Expr
    through: int
    derivation: str
    assumptions: s.logic.boolalg.Boolean = DOMAIN


def integer_polylog_example(rank=3, through=1):
    """Coefficient C_n=(1+e)/[(n+1)^(p-2)(n+1+e)]."""
    if rank<2:
        raise ValueError('rank must be at least 2')
    e,z=EPSILON,Z
    expression=s.hyper([1]*(rank-1)+[1+e],[2]*(rank-2)+[2+e],z)
    expected=s.series((1+e)*sum((-e)**k*s.polylog(rank-1+k,z)
                                for k in range(through+1))/z,e,0,through+1).removeO()
    return Example(f'integer-{rank}F{rank-1}',expression,expected,through,
                   'Cancel factorials in the Pochhammer coefficient; expand 1/(n+1+epsilon).')


def matched_half_example(through=1):
    """Coefficient C_n=(1+2e)/[(2n+1)(2n+1+2e)]."""
    e,z=EPSILON,Z
    half=s.Rational(1,2)
    expression=s.hyper([1,half,half+e],[half+1,half+1+e],z)
    odd=lambda weight:(s.polylog(weight,s.sqrt(z))-s.polylog(weight,-s.sqrt(z)))/(2*s.sqrt(z))
    expected=s.series((1+2*e)*sum((-2*e)**k*odd(k+2) for k in range(through+1)),
                      e,0,through+1).removeO()
    return Example('matched-half-3F2',expression,expected,through,
                   'Cancel half-Pochhammer ratios; project the polylog series onto odd indices.')


def singular_endpoint_example(through=1):
    """Convergent for e<0: 2F1(1,1;2-e;1)=(e-1)/e."""
    e=EPSILON
    expression=s.gamma(e)*s.hyper([1,1],[2-e],1)
    expected=s.series(s.gamma(e)*(e-1)/e,e,0,through+1).removeO()
    return Example('Gauss-endpoint-times-Gamma',expression,expected,through,
                   'Gauss summation for -epsilon>0 plus Gamma recurrence; combined double pole.',
                   e<0)


def examples():
    return [integer_polylog_example(),matched_half_example(),singular_endpoint_example()]


def classical_reduction(expression):
    """Only shuffle log extraction and H(0^k,-1;x)=-Li_(k+1)(-x).

    The latter follows by differentiating and matching the zero-origin value.
    It is an exact common-domain identity, not a numerical fit.
    """
    result=extract_hpl_logs(expression)
    replacements={h:-s.polylog(len(h.word),-h.argument) for h in result.atoms(HPL)
                  if h.word[-1]==-1 and all(letter==0 for letter in h.word[:-1])}
    return s.expand_log(result.xreplace(replacements))


def run():
    for case in examples():
        result=expand_epsilon(case.expression,EPSILON,through=case.through,
                              assumptions=case.assumptions,config=NativeConfig(use_cache=False))
        assert s.simplify(s.expand(classical_reduction(result.polynomial)-case.expected))==0
        print(case.name, ': exact comparison passed')
        print('Native:',result.polynomial)
        print('Independent formula:',case.expected)
        print('Derivation:',case.derivation)


if __name__=='__main__':
    run()
