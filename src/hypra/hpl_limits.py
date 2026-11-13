"""Combined finite z->1- limits using exact HPL reflection or power/log series.

This is not termwise endpoint evaluation: rational/logarithmic singularities
are retained until all coefficients at the reflected origin have combined.
"""
import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import _as_hpl, _expansion_size, _origin_constant
from .hpl_transform import reflect_hpl
from .symbolic import HPL


def hpl_limit_at_one(expression, z):
    """Prove the combined finite real endpoint, or raise on unsupported data.

    Coefficients must be rational in z, with ordinary HPLs/logs/polylogs at z.
    Reflection supports binary words of weight<=8. Ternary words use exact
    near-one series (individual weight<=8, total polynomial weight<=12,
    rational pole order<=16, numeric rational-function coefficients). No complex branch choice or
    divergent regularization is inferred. A failure is not an infinite value.
    """
    if not isinstance(expression,s.Expr) or not isinstance(z,s.Symbol):
        raise TypeError('Expected a SymPy expression and symbolic argument')
    if expression.has(s.Float,s.nan,s.zoo,s.oo,-s.oo):
        raise UnsupportedExpressionError('Endpoint expression must be exact and finite')
    transformed=_as_hpl(expression,z)
    if any(-1 in h.word for h in transformed.atoms(HPL)):
        return _ternary_limit(transformed, z)
    t=s.Dummy('endpoint_distance',positive=True)
    replacements={h:reflect_hpl(h.word,t,domain='0<z<1') for h in transformed.atoms(HPL)}
    reflected=transformed.xreplace(replacements).subs(z,1-t)
    return _origin_constant(reflected,t)


def _ternary_limit(expression, z):
    """Combine exact near-one Laurent/log coefficients before testing finiteness."""
    from .hpl_series_one import hpl_series_at_one

    if _expansion_size(expression) > 4096:
        raise UnsupportedExpressionError('Endpoint expression exceeds the expansion budget')
    atoms = sorted(expression.atoms(HPL), key=s.default_sort_key)
    generators = [s.Dummy('endpoint_hpl') for _ in atoms]
    replaced = expression.xreplace(dict(zip(atoms, generators, strict=True)))
    try:
        polynomial = s.Poly(replaced, *generators)
    except s.PolynomialError as exc:
        raise UnsupportedExpressionError('Ternary endpoint limit requires polynomial HPL dependence') from exc
    t, T = s.Dummy('endpoint_distance', positive=True), s.Dummy('endpoint_log')
    terms, pole = [], 0
    for powers, coefficient in polynomial.terms():
        if sum(p*len(h.word) for p,h in zip(powers,atoms,strict=True)) > 12:
            raise UnsupportedExpressionError('Endpoint expression exceeds total HPL weight twelve')
        if coefficient.free_symbols - {z} or not coefficient.is_rational_function(z):
            raise UnsupportedExpressionError('Ternary endpoint coefficients must be numeric rational functions of z')
        coefficient = s.cancel(coefficient.subs(z, 1-t))
        numerator, denominator = coefficient.as_numer_denom()
        try:
            num, den = s.Poly(numerator,t), s.Poly(denominator,t)
        except s.PolynomialError as exc:
            raise UnsupportedExpressionError('Unresolved rational endpoint coefficient') from exc
        if num.is_zero:
            continue
        nmin, dmin = min(m[0] for m in num.monoms()), min(m[0] for m in den.monoms())
        if den.nth(dmin).is_zero is not False:
            raise UnsupportedExpressionError('Could not prove the endpoint denominator nonzero')
        pole = max(pole, dmin-nmin)
        terms.append((powers, coefficient))
    if pole > 16:
        raise UnsupportedExpressionError('Endpoint rational pole order exceeds sixteen')
    expansions = [hpl_series_at_one(h.word,t,through=pole).removeO().subs(s.log(t),-T)
                  for h in atoms]
    result = s.S.Zero
    for powers, coefficient in terms:
        result += s.series(coefficient,t,0,1).removeO()*s.prod(
            value**power for value,power in zip(expansions,powers,strict=True))
    _expansion_size(result)
    # Multiply by a sufficient t power to inspect an ordinary polynomial.
    combined = s.Poly(s.expand(result*t**pole), t, T)
    constant = s.S.Zero
    for (degree, logarithms), coefficient in combined.terms():
        if degree > pole:
            continue
        coefficient = s.simplify(coefficient)
        if degree < pole or logarithms:
            if coefficient != 0:
                raise UnsupportedExpressionError('Endpoint diverges or cancellation is not proved')
        else:
            constant += coefficient
    return s.simplify(constant)
