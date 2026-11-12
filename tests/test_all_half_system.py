"""All-half 3F2 checks against the defining coefficient product, not an oracle."""
import pytest
import sympy as s

from hypra.half_system import canonical_half_series
from hypra.hpl_algebra import hpl_series

z,e = s.symbols('z e')
t = s.Symbol('t', positive=True)
half = s.Rational(1,2)


@pytest.mark.parametrize('alpha,beta', [([1,2,-1],[3,-2]), ([0,1,2],[1,0]), ([1,2,-1,0],[3,-2,1])])
def test_all_half_exact_pochhammer_coefficients(alpha,beta):
    order, degree = 2, 5
    result = canonical_half_series([half]*len(alpha),alpha,[half]*len(beta),beta,z,e,order)
    # Remove the known exact gauge before the rational change of variable.
    y = s.expand(result*s.sqrt(1-z)).xreplace({s.I*s.sqrt(z/(1-z)):t})
    assert not y.has(z)
    actual = s.series(hpl_series(y,t,degree)*s.sqrt(1-t*t),t,0,degree).removeO().expand()
    zt=t*t/(t*t-1)
    expected=s.S.Zero
    for n in range(3):
        c=s.prod(s.rf(half+a*e,n) for a in alpha)/(s.factorial(n)*s.prod(s.rf(half+b*e,n) for b in beta))
        expected += s.series(c,e,0,order+1).removeO()*s.series(zt**n,t,0,degree).removeO()
    assert s.expand(actual-expected)==0


def test_all_half_zeroth_order_and_exact_cancellation():
    assert canonical_half_series([half]*3,[1,2,3],[half]*2,[4,5],z,e,0)==1/s.sqrt(1-z)
    # Equal pairs cancel identically to 1F0(1/2+e;;z).
    result=canonical_half_series([half]*3,[1,2,3],[half]*2,[2,3],z,e,1)
    from hypra.symbolic import reduce_known_hpl
    normalized=reduce_known_hpl(s.expand(result*s.sqrt(1-z)))
    # HPL continuation may intentionally retain split complex logs; check
    # the exact Taylor germ rather than force a branch-changing log identity.
    yt=normalized.xreplace({s.I*s.sqrt(z/(1-z)):t})
    if yt.has(z):
        # Use unreduced words for the independent Taylor check.
        yt=s.expand(result*s.sqrt(1-z)).xreplace({s.I*s.sqrt(z/(1-z)):t})
    assert s.expand(hpl_series(yt,t,5)-(1-e*t**2-e*t**4/2))==0
