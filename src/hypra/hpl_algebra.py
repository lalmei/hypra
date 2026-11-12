"""Exact shuffle algebra and definite rational-HPL integration.

The domain is the real interval 0 < z < 1, with the HPL conventions in
``hypra.symbolic``. Rational poles may occur only at 0, 1, and -1. An origin
limit is computed for the complete primitive, so cancellations between
separately divergent summands are retained. No unevaluated Integral or guessed
integration constant is returned.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

import sympy as s

from .errors import UnsupportedExpressionError
from .symbolic import HPL

_MAX_WEIGHT = 24
_MAX_TERMS = 16384
_MAX_POLE = 32
_LOG = s.Dummy("log_z", real=True)


def _h(word: tuple[int, ...], z: s.Symbol) -> s.Expr:
    return HPL(s.Tuple(*word), z) if word else s.S.One


def _word(word) -> tuple[int, ...]:
    if not isinstance(word, (tuple, list, s.Tuple)):
        raise UnsupportedExpressionError("An HPL word must be a sequence of -1, 0, 1")
    if any(
        isinstance(x, bool) or x not in (-1, 0, 1) or not isinstance(x, (int, s.Integer))
        for x in word
    ):
        raise UnsupportedExpressionError("An HPL word must contain exact letters -1, 0, 1")
    result = tuple(map(int, word))
    if len(result) > _MAX_WEIGHT:
        raise UnsupportedExpressionError(f"HPL integration supports weight at most {_MAX_WEIGHT}")
    return result


@lru_cache(maxsize=1024)
def _shuffle(left: tuple[int, ...], right: tuple[int, ...]):
    if len(left) + len(right) > _MAX_WEIGHT:
        raise UnsupportedExpressionError(f"Shuffle weight exceeds {_MAX_WEIGHT}")
    if not left:
        return ((right, 1),)
    if not right:
        return ((left, 1),)
    counts = defaultdict(int)
    for word, count in _shuffle(left[1:], right):
        counts[(left[0], *word)] += count
    for word, count in _shuffle(left, right[1:]):
        counts[(right[0], *word)] += count
    if len(counts) > _MAX_TERMS:
        raise UnsupportedExpressionError("Shuffle expression exceeds the supported term budget")
    return tuple(sorted(counts.items()))


def shuffle_words(left, right) -> dict[tuple[int, ...], int]:
    """Return the shuffle words with their exact integer multiplicities."""
    return dict(_shuffle(_word(left), _word(right)))


def shuffle_product(left: HPL, right: HPL) -> s.Expr:
    """Multiply two same-argument HPLs using the shuffle identity."""
    if not isinstance(left, HPL) or not isinstance(right, HPL):
        raise TypeError("shuffle_product requires two HPL objects")
    if left.argument != right.argument:
        raise UnsupportedExpressionError("Shuffle products require a common HPL argument")
    return s.Add(
        *(
            count * _h(word, left.argument)
            for word, count in _shuffle(_word(left.word), _word(right.word))
        )
    )


def _input(expression, z):
    if not isinstance(z, s.Symbol):
        raise TypeError("The integration variable must be a SymPy Symbol")
    if isinstance(expression, bool):
        raise TypeError("A boolean is not an integration expression")
    if isinstance(expression, int):
        expression = s.Integer(expression)
    if not isinstance(expression, s.Expr):
        raise TypeError("The integrand must be an exact SymPy expression, not source text")
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("The integrand must be exact and finite")
    return expression


def _as_hpl(expression, z):
    """Convert only identities with matching principal real branch conventions."""
    replacements = {}
    for atom in expression.atoms(s.log):
        if not atom.has(z):
            continue
        argument = atom.args[0]
        if argument == z:
            replacements[atom] = _h((0,), z)
        elif argument == 1 - z:
            replacements[atom] = -_h((1,), z)
        elif argument == 1 + z:
            replacements[atom] = _h((-1,), z)
        else:
            raise UnsupportedExpressionError("Only log(z), log(1-z), log(1+z) are integrated")
    for atom in expression.atoms(s.polylog):
        if not atom.has(z):
            continue
        order, argument = atom.args
        if not isinstance(order, s.Integer) or not 1 <= order <= _MAX_WEIGHT:
            raise UnsupportedExpressionError(
                "Only positive-integer classical polylog orders are supported"
            )
        if argument == z:
            replacements[atom] = _h((0,) * (int(order) - 1) + (1,), z)
        elif argument == -z:
            replacements[atom] = -_h((0,) * (int(order) - 1) + (-1,), z)
        else:
            raise UnsupportedExpressionError("Only Li_n(z) and Li_n(-z) are converted to HPLs")
    for atom in expression.atoms(HPL):
        _word(atom.word)
        if atom.argument != z:
            raise UnsupportedExpressionError(
                "HPL integration requires the integration variable as argument"
            )
    return expression.xreplace(replacements)


def _expansion_size(expression):
    """Conservative pre-expansion budget, avoiding a giant multinomial first."""
    if isinstance(expression, HPL):
        return 1
    if expression.is_Add:
        size = sum(_expansion_size(arg) for arg in expression.args)
    elif expression.is_Mul:
        size = 1
        for arg in expression.args:
            size *= _expansion_size(arg)
            if size > _MAX_TERMS:
                break
    elif expression.is_Pow and isinstance(expression.exp, s.Integer) and expression.exp > 0:
        if expression.exp > 128:
            raise UnsupportedExpressionError("Polynomial power exceeds the supported budget")
        base_size = _expansion_size(expression.base)
        size = base_size ** int(expression.exp)
    else:
        size = 1
    if size > _MAX_TERMS:
        raise UnsupportedExpressionError("Expression expansion exceeds the supported term budget")
    return size


def _linear_terms(expression, z):
    """Map words to rational coefficients, after exact shuffle expansion."""
    expression = _as_hpl(_input(expression, z), z)
    _expansion_size(expression)
    expression = s.expand(expression, log=False, power_base=False, power_exp=False)
    if len(s.Add.make_args(expression)) > _MAX_TERMS:
        raise UnsupportedExpressionError("HPL expression exceeds the supported term budget")
    result = defaultdict(lambda: s.S.Zero)
    for term in s.Add.make_args(expression):
        coefficient = s.S.One
        words = {(): s.S.One}
        for factor in s.Mul.make_args(term):
            base, exponent = factor.as_base_exp()
            if isinstance(base, HPL) and base.argument == z:
                if not isinstance(exponent, s.Integer) or exponent < 0:
                    raise UnsupportedExpressionError(
                        "Only nonnegative integer powers of HPLs are supported"
                    )
                word = _word(base.word)
                if len(word) * exponent > _MAX_WEIGHT:
                    raise UnsupportedExpressionError("HPL product exceeds supported total weight")
                for _ in range(int(exponent)):
                    next_words = defaultdict(lambda: s.S.Zero)
                    for existing, multiplier in words.items():
                        for combined, count in _shuffle(existing, word):
                            next_words[combined] += multiplier * count
                    words = next_words
                    if len(words) > _MAX_TERMS:
                        raise UnsupportedExpressionError(
                            "HPL product exceeds the supported term budget"
                        )
            else:
                if any(atom.argument.has(z) for atom in factor.atoms(HPL)):
                    raise UnsupportedExpressionError(
                        "HPLs must occur polynomially, not inside other functions"
                    )
                coefficient *= factor
        coefficient = s.cancel(coefficient)
        if coefficient.is_rational_function(z) is not True:
            raise UnsupportedExpressionError(
                "HPL coefficients must be rational in the integration variable"
            )
        for word, multiplier in words.items():
            result[word] += coefficient * multiplier
    return {word: s.cancel(coefficient) for word, coefficient in result.items() if coefficient != 0}


def shuffle_expand(expression, z: s.Symbol) -> s.Expr:
    """Linearize HPL products, powers and supported logs/polylogs exactly.

    The result is a sum of rational coefficients times single HPLs. Matching
    logarithms use the real principal branch on 0 < z < 1. No argument
    transformations, numeric fitting or branch-changing log identities occur.
    """
    return s.Add(
        *(coefficient * _h(word, z) for word, coefficient in _linear_terms(expression, z).items())
    )


def _pole_multiplicities(denominator, z):
    poly = s.Poly(denominator, z, domain="EX")
    multiplicities = {}
    for pole in (0, 1, -1):
        factor = s.Poly(z - pole, z, domain="EX")
        count = 0
        while poly.degree() > 0:
            quotient, remainder = s.div(poly, factor)
            if not remainder.is_zero:
                break
            poly = quotient
            count += 1
            if count > _MAX_POLE:
                raise UnsupportedExpressionError(f"Rational pole order exceeds {_MAX_POLE}")
        multiplicities[pole] = count
    if poly.degree() != 0:
        raise UnsupportedExpressionError(
            "Rational HPL integration supports poles only at 0, 1 and -1"
        )
    return multiplicities


@lru_cache(maxsize=512)
def _split_rational(rational, z):
    """Return R0 and kernel residues such that R=R0'+sum c_a*f_a."""
    rational = s.cancel(rational)
    if rational == 0:
        return s.S.Zero, ()
    numerator, denominator = s.fraction(rational)
    multiplicities = _pole_multiplicities(denominator, z)
    polynomial, _ = s.div(s.Poly(numerator, z, domain="EX"), s.Poly(denominator, z, domain="EX"))
    primitive = s.Add(
        *(
            coefficient * z ** (power[0] + 1) / (power[0] + 1)
            for power, coefficient in polynomial.terms()
        )
    )
    residues = {}
    for pole, multiplicity in multiplicities.items():
        if not multiplicity:
            continue
        regular = s.cancel(rational * (z - pole) ** multiplicity)
        for order in range(1, multiplicity + 1):
            derivative_order = multiplicity - order
            coefficient = s.cancel(s.diff(regular, z, derivative_order)).subs(
                z, pole
            ) / s.factorial(derivative_order)
            if coefficient.has(s.nan, s.zoo, s.oo, -s.oo, z):
                raise UnsupportedExpressionError(
                    "Could not establish an exact rational partial fraction"
                )
            if order == 1:
                residues[pole] = -coefficient if pole == 1 else coefficient
            else:
                primitive += coefficient * (z - pole) ** (1 - order) / (1 - order)
    return s.cancel(primitive), tuple(
        (pole, coefficient) for pole, coefficient in residues.items() if coefficient != 0
    )


@lru_cache(maxsize=1024)
def _primitive(rational, word, z):
    rational_primitive, residues = _split_rational(s.cancel(rational), z)
    result = rational_primitive * _h(word, z)
    result += s.Add(*(coefficient * _h((pole, *word), z) for pole, coefficient in residues))
    if word and rational_primitive != 0:
        kernel = {0: 1 / z, 1: 1 / (1 - z), -1: 1 / (1 + z)}[word[0]]
        result -= _primitive(s.cancel(rational_primitive * kernel), word[1:], z)
    return result


def _log_primitive(polynomial, power):
    """Polynomial Q with d(z^power Q(log z))/dz=z^(power-1) P(log z)."""
    if power <= 0:
        raise UnsupportedExpressionError("Internal HPL local-series power is not positive")
    result = s.S.Zero
    for (degree,), coefficient in s.Poly(polynomial, _LOG).terms():
        for derivative in range(degree + 1):
            result += (
                coefficient
                * (-1) ** derivative
                * s.factorial(degree)
                / s.factorial(degree - derivative)
                * _LOG ** (degree - derivative)
                / power ** (derivative + 1)
            )
    return s.expand(result)


@lru_cache(maxsize=1024)
def _local_hpl(word, through):
    """Finite exact expansion in z, with polynomial coefficients in log(z)."""
    if not word:
        return ((0, s.S.One),)
    if all(letter == 0 for letter in word):
        return ((0, _LOG ** len(word) / s.factorial(len(word))),)
    if through <= 0:
        return ()
    head, tail = word[0], word[1:]
    result = defaultdict(lambda: s.S.Zero)
    for power, polynomial in _local_hpl(tail, through if head == 0 else through - 1):
        if head == 0:
            if power <= 0:
                raise UnsupportedExpressionError("Unresolved logarithmic origin term")
            result[power] += _log_primitive(polynomial, power)
        else:
            for kernel_power in range(through - power):
                target_power = power + kernel_power + 1
                sign = (-1) ** kernel_power if head == -1 else 1
                result[target_power] += sign * _log_primitive(polynomial, target_power)
    return tuple((power, s.expand(polynomial)) for power, polynomial in sorted(result.items()))


def _origin_constant(expression, z):
    """Prove a finite origin limit after combining all Laurent-log terms."""
    powers = defaultdict(lambda: s.S.Zero)
    for word, rational in _linear_terms(expression, z).items():
        if rational == 0:
            continue
        numerator, denominator = s.fraction(s.cancel(rational))
        del numerator
        pole_order = _pole_multiplicities(denominator, z)[0]
        regular = s.cancel(rational * z**pole_order)
        hpl_series = _local_hpl(word, pole_order)
        for j in range(pole_order + 1):
            rational_coefficient = s.cancel(s.diff(regular, z, j)).subs(z, 0) / s.factorial(j)
            rational_power = j - pole_order
            for hpl_power, polynomial in hpl_series:
                power = rational_power + hpl_power
                if power <= 0:
                    powers[power] += rational_coefficient * polynomial
    origin = s.S.Zero
    for power, polynomial in powers.items():
        polynomial = s.Poly(s.expand(polynomial), _LOG)
        for (log_power,), coefficient in polynomial.terms():
            coefficient = s.simplify(coefficient)
            if coefficient == 0:
                continue
            if power < 0 or log_power > 0:
                raise UnsupportedExpressionError(
                    "The origin integral diverges or its convergence cannot be proved"
                )
            origin += coefficient
    return s.simplify(origin)


def primitive_hpl(expression, z: s.Symbol) -> s.Expr:
    """An exact local primitive of rational-HPL products, up to a constant.

    Rational poles must be at 0,±1. Unlike integrate_from_zero, this does not
    require a finite origin value and does not assign a divergent endpoint.
    Principal branches on real 0<z<1 are retained; no definite integral or
    contour continuation is inferred from this local antiderivative.
    """
    linear = _linear_terms(expression, z)
    primitive = s.Add(*(_primitive(rational, word, z) for word, rational in linear.items()))
    return shuffle_expand(primitive, z)


def integrate_from_zero(expression, z: s.Symbol) -> s.Expr:
    """Compute integral_0^z R(t) HPL(word,t) dt, including shuffle products.

    Rational coefficients may have poles only at 0, 1, -1; polynomial parts,
    repeated poles, exact parameter coefficients and finite cancellations at
    zero are supported. The result applies for 0 < z < 1 with real principal
    branches. Unsupported structure, excessive resources or an unproved finite
    origin limit raises UnsupportedExpressionError; no divergent terms are
    silently regularized. The algorithm recursively lowers HPL word weight.
    """
    primitive = primitive_hpl(expression, z)
    origin = _origin_constant(primitive, z)
    return primitive - origin


def inverse_theta_plus(expression, z: s.Symbol, n: int) -> s.Expr:
    """Regular inverse of theta+n: z^-n integral_0^z t^(n-1) f(t) dt.

    Here theta=z*d/dz and n is a positive integer. The vanishing integration
    constant selects the solution without an added homogeneous z^-n term.
    """
    if isinstance(n, bool) or not isinstance(n, (int, s.Integer)) or not 1 <= n <= _MAX_POLE:
        raise UnsupportedExpressionError(
            f"inverse_theta_plus requires an integer 1 <= n <= {_MAX_POLE}"
        )
    return integrate_from_zero(z ** (n - 1) * _input(expression, z), z) / z**n


def hpl_series(expression, z: s.Symbol, order: int) -> s.Expr:
    """Return exact Laurent/log terms of argument degree strictly below order.

    This is an independent local expansion at z=0 on the positive real side,
    obtained by integrating the defining HPL kernels coefficientwise. It
    includes negative rational powers and powers of log(z). The return value
    is a finite truncation, not an equality to the full function and not a
    falsely asserted O(z**order) remainder (which may carry logarithms).
    """
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)) or not 1 <= order <= 32:
        raise UnsupportedExpressionError("hpl_series requires an integer 1 <= order <= 32")
    powers = defaultdict(lambda: s.S.Zero)
    for word, rational in _linear_terms(expression, z).items():
        if rational == 0:
            continue
        _, denominator = s.fraction(s.cancel(rational))
        pole_order = _pole_multiplicities(denominator, z)[0]
        through = int(order) - 1 + pole_order
        if through > 48:
            raise UnsupportedExpressionError("Local HPL series exceeds the supported power budget")
        regular = s.cancel(rational * z**pole_order)
        local = _local_hpl(word, through)
        for j in range(through + 1):
            coefficient = s.cancel(s.diff(regular, z, j)).subs(z, 0) / s.factorial(j)
            rational_power = j - pole_order
            for hpl_power, polynomial in local:
                power = rational_power + hpl_power
                if power < order:
                    powers[power] += coefficient * polynomial
    return s.Add(
        *(
            z**power * s.expand(polynomial).subs(_LOG, s.log(z))
            for power, polynomial in sorted(powers.items())
        )
    )
