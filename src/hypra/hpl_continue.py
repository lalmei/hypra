"""Explicit radial numerical HPL continuation by bounded Taylor ODE steps.

The mathematical truncation bound excludes floating-point roundoff. A second
higher-precision run is a stability check, not an interval certificate.
"""
from dataclasses import dataclass

import mpmath as mp
import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_numeric import (
    _DECIMAL,
    HPLNumericalConvergenceError,
    _word,
    evaluate_hpl,
)


@dataclass(frozen=True)
class HPLContinuationResult:
    value: object
    tail_bound: object
    stability_gap: object
    steps: int
    max_degree: int
    dps: int
    working_dps: int
    method: str = 'triangular-suffix-ODE-Taylor'
    domain: str = 'principal branch off real cuts; straight radial path from origin germ'


def _target(value, ctx):
    if isinstance(value, bool):
        raise TypeError('A boolean is not a numerical argument')
    if isinstance(value, s.Expr):
        if value.free_symbols or value.is_finite is not True:
            raise UnsupportedExpressionError('Use a finite numerical argument')
        real, imag = value.as_real_imag()
        z = ctx.mpc(str(s.N(real, ctx.dps+10)), str(s.N(imag, ctx.dps+10)))
    elif isinstance(value, str):
        if len(value)>4096 or not _DECIMAL.fullmatch(value):
            raise ValueError('Use a decimal string, not executable text')
        z = ctx.mpc(value)
    else:
        try:
            z = ctx.mpc(value)
        except (TypeError, ValueError) as exc:
            raise TypeError('Use an exact SymPy number or a numerical scalar') from exc
    if not ctx.isfinite(z) or z == 0 or (ctx.im(z)==0 and not 0<ctx.re(z)<1):
        raise UnsupportedExpressionError('Radial continuation excludes zero, endpoints and real cuts (-infinity,0] and [1,infinity)')
    return z


def _run(word, argument, dps, working, max_steps, max_degree, max_terms):
    ctx=mp.mp.clone()
    ctx.dps=working
    target=_target(argument,ctx)
    c=target*min(ctx.one,ctx.mpf('0.25')/abs(target))
    # Exact decimal data avoids cross-context mpc type/coercion limitations.
    initial=s.Rational(str(ctx.re(c)))+s.I*s.Rational(str(ctx.im(c)))
    values=[]
    error=ctx.zero
    for j in range(len(word)):
        evaluated=evaluate_hpl(word[j:],initial,dps=working,max_terms=max_terms,domain='principal-disk')
        values.append(ctx.mpc(evaluated.value))
        error=max(error,ctx.mpf(evaluated.tail_bound))
    values.append(ctx.one)
    steps=degree_used=0
    local_tolerance=ctx.power(10,-dps-15)/max_steps
    while c != target:
        if steps >= max_steps:
            raise HPLNumericalConvergenceError('Radial continuation exceeded max_steps')
        distance=min(abs(c),abs(c-1),abs(c+1))
        h=target-c
        if abs(h)>distance/4:
            h=h*(distance/(4*abs(h)))
        radius=distance/2
        ratio=abs(h)/radius
        majorant=max(map(abs,values))*ctx.exp(radius/(distance-radius))
        degree=max(1,int(ctx.ceil(ctx.log(local_tolerance*(1-ratio)/majorant)/ctx.log(ratio)))-1)
        if degree>max_degree:
            raise HPLNumericalConvergenceError('Radial continuation exceeded max_degree')
        coefficients=[[v]+[ctx.zero]*degree for v in values]
        for j in range(len(word)-1,-1,-1):
            letter=word[j]
            denominator=c if letter==0 else (1-c if letter==1 else 1+c)
            sign=1 if letter==1 else -1
            kernels=[ctx.one/denominator]
            for _ in range(degree-1):
                kernels.append(kernels[-1]*sign/denominator)
            for n in range(degree):
                coefficients[j][n+1]=sum(kernels[k]*coefficients[j+1][n-k] for k in range(n+1))/(n+1)
        values=[ctx.polyval(list(reversed(row)),h) for row in coefficients]
        local_error=majorant*ratio**(degree+1)/(1-ratio)
        error=error*ctx.exp(abs(h)/(distance-abs(h)))+local_error
        c=target if h==target-c else c+h
        steps+=1
        degree_used=max(degree_used,degree)
    if error>ctx.power(10,-dps)*max(ctx.one,abs(values[0])):
        raise HPLNumericalConvergenceError('Propagated truncation bound exceeds requested tolerance')
    return values[0],error,steps,degree_used


def continue_hpl(word, argument, *, domain, path, dps=50, max_steps=4096,
                 max_degree=1024, max_terms=100000):
    """Numerically continue along an explicitly requested principal radial path.

    Required choices: domain='principal-cut-plane', path='radial'. Real input
    is permitted only between zero and one. Complex paths must be off the real
    cuts. No boundary value, winding path or exact symbolic identity is inferred.
    """
    if domain!='principal-cut-plane' or path!='radial':
        raise UnsupportedExpressionError('Specify domain="principal-cut-plane", path="radial" explicitly')
    word=_word(word)
    for name,value,upper in [('dps',dps,500),('max_steps',max_steps,100000),
                             ('max_degree',max_degree,4096),('max_terms',max_terms,1000000)]:
        if isinstance(value,bool) or not isinstance(value,int) or not 1<=value<=upper:
            raise ValueError(f'{name} must be an integer between 1 and {upper}')
    first=_run(word,argument,dps,dps+25,max_steps,max_degree,max_terms)
    second=_run(word,argument,dps,dps+45,max_steps,max_degree,max_terms)
    ctx=mp.mp.clone()
    ctx.dps=dps+45
    gap=abs(ctx.mpc(first[0])-ctx.mpc(second[0]))
    if gap>ctx.power(10,-dps)*max(ctx.one,abs(second[0])):
        raise HPLNumericalConvergenceError('Higher-precision continuation did not stabilize')
    return HPLContinuationResult(second[0],second[1],gap,second[2],second[3],dps,dps+45)
