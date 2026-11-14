"""Power/log limits of the rationalized hypergeometric coefficients at one."""

from functools import lru_cache

import sympy as s

from ..errors import UnsupportedExpressionError
from ..hpl.algebra import _LOG, _expansion_size, _local_hpl
from ..hpl.inversion import invert_hpl
from ..hpl.limits import hpl_limit_at_one
from ..hpl.series_one import hpl_series_at_one
from ..symbolic import HPL


@lru_cache(maxsize=256)
def _inverted(word, variable):
    return invert_hpl(word, variable, boundary="upper", domain="0<z<1")


def unit_coefficient_limit(expression, z):
    """Combine finite Laurent/log terms using the positive sqrt(1-z) coordinate."""
    if all(atom.argument == z for atom in expression.atoms(HPL)) and all(
        atom.args[1] in (z, -z) for atom in expression.atoms(s.polylog)
    ) and not expression.has(s.sqrt(1-z), s.sqrt(z)):
        return hpl_limit_at_one(expression, z)

    t = s.Dummy("endpoint_root", positive=True)
    y = s.Dummy("inverse_argument", positive=True)
    coordinates = (
        z, s.sqrt(z), (1-s.sqrt(1-z))/(1+s.sqrt(1-z)), s.I*s.sqrt(z/(1-z)),
    )
    # Restore the HPL form of the simple polylogs emitted by the kernels.
    polylogs = {}
    for atom in expression.atoms(s.polylog):
        weight, argument = atom.args
        for coordinate in coordinates:
            if s.cancel(argument-coordinate) == 0:
                polylogs[atom] = HPL((0,)*(int(weight)-1)+(1,), coordinate)
                break
            if s.cancel(argument+coordinate) == 0:
                polylogs[atom] = -HPL((0,)*(int(weight)-1)+(-1,), coordinate)
                break
        else:
            raise UnsupportedExpressionError("Unrecognized polylogarithm argument at z=1")
    expression = expression.xreplace(polylogs)
    atoms = sorted(expression.atoms(HPL), key=s.default_sort_key)
    generators = [s.Dummy("endpoint_word") for _ in atoms]
    _expansion_size(expression)
    replaced = expression.xreplace(dict(zip(atoms, generators, strict=True)))
    try:
        terms = s.Poly(replaced, *generators).terms() if generators else [((), replaced)]
    except s.PolynomialError as exc:
        raise UnsupportedExpressionError("Endpoint requires polynomial HPL dependence") from exc

    coefficients = []
    depth = 0
    for powers, coefficient in terms:
        coefficient = coefficient.subs(z, 1-t**2)
        leading = coefficient.as_leading_term(t)
        _, power = leading.as_coeff_exponent(t)
        if leading == 0:
            continue
        if not power.is_Integer or power < -24:
            raise UnsupportedExpressionError("Unsupported rationalized endpoint power")
        depth = max(depth, -int(power))
        if power <= 0:
            coefficients.append((powers, coefficient, int(power)))

    local = []
    for atom in atoms:
        argument = s.simplify(atom.argument.subs(z, 1-t**2))
        if s.limit(argument, t, 0, dir="+") == 1:
            distance = s.Symbol("_unit_distance", positive=True)
            germ = hpl_series_at_one(atom.word, distance, through=depth).removeO()
            germ = germ.subs(distance, s.factor(1-argument))
        elif s.limit(argument/t**-1, t, 0, dir="+") == s.I:
            # The upper real inversion germ continues to the upper
            # imaginary axis, with inverse coordinate approaching from below.
            inverted = _inverted(tuple(atom.word), y)
            substitutions = {}
            for hpl in inverted.atoms(HPL):
                substitutions[hpl] = s.Add(*(
                    y**power * polynomial.subs(_LOG, s.log(y))
                    for power, polynomial in _local_hpl(tuple(hpl.word), depth)
                ))
            germ = inverted.xreplace(substitutions).subs(y, 1/argument)
        else:
            raise UnsupportedExpressionError("Unrecognized HPL coordinate at z=1")
        local.append(s.series(germ, t, 0, depth+1).removeO())

    result = s.S.Zero
    for powers, coefficient, valuation in coefficients:
        rational = s.series(coefficient, t, 0, 1).removeO().subs(s.log(t), _LOG)
        product = s.Poly(1, t)
        for germ, power in zip(local, powers, strict=True):
            if power:
                product *= s.Poly(germ.subs(s.log(t), _LOG), t)**power
                product = s.Poly.from_dict(
                    {degree: term for degree, term in product.terms() if degree[0] <= -valuation}, t
                )
        for (power,), term in s.Poly(s.expand(t**depth*rational), t).terms():
            if power <= depth:
                result += s.Add(*(term*product.nth(j)*t**(power+j-depth)
                                  for j in range(depth-power+1)))
    value = s.simplify(s.expand(result))
    if value.has(t, _LOG, s.oo, -s.oo, s.zoo, s.nan, s.Limit):
        raise UnsupportedExpressionError("Could not establish a finite rationalized endpoint")
    return value
