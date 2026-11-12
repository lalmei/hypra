"""Independent exact coefficient checks against defining Pochhammer series."""
import functools

import pytest
import sympy as s

from hypra.half_system import canonical_half_series
from hypra.symbolic import HPL

z, e = s.symbols('z e')
t = s.Symbol('t', positive=True)
half = s.Rational(1, 2)


@functools.lru_cache(None)
def hpl_taylor(word, degree):
    if not word:
        return s.S.One
    if all(letter == 0 for letter in word):
        return s.log(t)**len(word) / s.factorial(len(word))
    kernel = {0: 1/t, 1: 1/(1-t), -1: 1/(1+t)}[word[0]]
    integrand = s.series(kernel*hpl_taylor(word[1:], degree), t, 0, degree).removeO()
    return s.integrate(integrand, t)


def transformed(result, kind):
    if kind == 'A':
        argument = s.I*s.sqrt(z/(1-z))
        zt = t**2/(t**2-1)
    elif kind == 'B':
        argument = (1-s.sqrt(1-z))/(1+s.sqrt(1-z))
        zt = 4*t/(1+t)**2
    else:
        argument, zt = s.sqrt(z), t**2
    result = result.xreplace({argument: t})
    assert not result.has(z)
    return result, zt


@pytest.mark.parametrize('a,b,kind', [
    ([0, 0], [half], 'A'),
    ([0, 0, 0], [1, half], 'A'),
    ([0, half], [1], 'B'),
    ([half, 0, 0], [1, 1], 'B'),
    ([0, half], [half], 'C'),
])
@pytest.mark.parametrize('slopes', [(1, 2, 3), (-1, 0, 2)])
def test_exact_defining_series_coefficients(a, b, kind, slopes):
    alpha = [s.Integer(slopes[i % 3]) for i in range(len(a))]
    beta = [s.Integer(slopes[(i+2) % 3]) for i in range(len(b))]
    order, degree = 3, 5
    actual, zt = transformed(canonical_half_series(a, alpha, b, beta, z, e, order), kind)
    actual = actual.xreplace({h: hpl_taylor(tuple(h.word), degree) for h in actual.atoms(HPL)})
    actual = s.series(actual, t, 0, degree).removeO().expand()
    expected = s.S.Zero
    # n >= degree cannot affect t powers < degree for any of the three maps.
    for n in range(degree):
        coefficient = s.prod(s.rf(v+w*e, n) for v,w in zip(a,alpha)) / (
            s.factorial(n)*s.prod(s.rf(v+w*e, n) for v,w in zip(b,beta)))
        coefficient = s.series(coefficient, e, 0, order+1).removeO()
        expected += coefficient*s.series(zt**n, t, 0, degree).removeO()
    assert s.expand(actual-expected) == 0


def test_symbolic_slopes_and_zero_order():
    alpha,beta,gamma = s.symbols('alpha beta gamma')
    value = canonical_half_series([0, half], [alpha,beta], [half], [gamma], z,e,1)
    assert s.expand(value - 1 - e*alpha*(HPL(s.Tuple(1), s.sqrt(z))-HPL(s.Tuple(-1),s.sqrt(z)))) == 0
    assert canonical_half_series([0, half], [alpha,beta], [half], [gamma], z,e,0) == 1


def test_unsupported_patterns_and_validation():
    assert canonical_half_series([half], [1], [], [], z,e,2) is None
    assert canonical_half_series([half,half], [1,2], [1], [3], z,e,2) is None
    with pytest.raises(ValueError):
        canonical_half_series([0,half], [1], [1], [3], z,e,2)
    with pytest.raises(ValueError):
        canonical_half_series([0,half], [1,e], [1], [3], z,e,2)
    with pytest.raises(ValueError):
        canonical_half_series([0,half], [1,2], [1], [3], z,e,True)


def test_vanishing_operator_coefficient():
    value = canonical_half_series([0,0], [1,-1], [half], [0], z,e,3)
    assert s.expand(value).coeff(e, 1) == 0
    assert s.expand(value).coeff(e, 3) == 0
