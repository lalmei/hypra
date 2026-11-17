"""Named-function reductions generated from exact iterated-integral identities.

The weight-four/five Nielsen relations needed by HypExp can be obtained by
linear algebra, without a formula table. Expand weighted monomials of classical
polylogarithms and logarithms into a common binary HPL alphabet. Shuffle their
products, and solve for the desired Nielsen word. Constant terms from the
argument transformations are retained as separate rows, fixing integration
constants and the specified cut side along with the nonconstant terms.
"""

from functools import lru_cache
from typing import cast

import sympy as s

from ._typing import Boundary, Interval
from .errors import UnsupportedExpressionError
from .hpl.algebra import _linear_terms
from .nielsen import NielsenPolylog, reduce_nielsen
from .polylogs import transform_polylog_arguments

_Z = s.Symbol("polylog_reduction_argument", positive=True)
_PHASE, _ZETA3, _ZETA5 = s.symbols("polylog_phase polylog_zeta3 polylog_zeta5")
_ARGUMENTS = (
    _Z,
    1 - _Z,
    _Z / (_Z - 1),
    1 - 1 / _Z,
    1 / _Z,
    1 / (1 - _Z),
)
_GENERATORS = (
    (s.log(_Z), 1),
    (s.log(1 - _Z), 1),
    (s.I * s.pi, 1),
    (s.polylog(2, _Z), 2),
    (s.polylog(3, _Z), 3),
    (s.polylog(3, 1 - _Z), 3),
    (s.zeta(3), 3),
    *((s.polylog(weight, argument), weight) for weight in (4, 5) for argument in _ARGUMENTS[:3]),
    (NielsenPolylog(3, 2, _Z), 5),
    (s.zeta(5), 5),
)


@lru_cache(maxsize=128)
def _monomials(weight: int, offset: int = 0) -> tuple[s.Expr, ...]:
    if weight == 0:
        return (s.S.One,)
    if offset == len(_GENERATORS):
        return ()
    expression, degree = _GENERATORS[offset]
    return tuple(
        expression**power * tail
        for power in range(weight // degree + 1)
        for tail in _monomials(weight - power * degree, offset + 1)
    )


def _features(expression: s.Expr, boundary: Boundary = "upper"):
    """Map each HPL word and formal constant monomial to a rational value.

    Use i*pi rather than pi as the weight-one constant. The coefficients of
    principal cut phases are then rational too. Equating these formal rows is
    sufficient for an identity; no transcendence or constant-independence
    conjecture is needed. In particular, the empty HPL word is never dropped.
    """
    transformed = transform_polylog_arguments(expression, _Z, domain="0<z<1", boundary=boundary)
    features = {}
    for word, coefficient in _linear_terms(transformed, _Z).items():
        if coefficient.has(_Z):
            raise UnsupportedExpressionError("Named-function basis has a nonconstant coefficient")
        formal = coefficient.xreplace({s.pi: -s.I * _PHASE, s.zeta(3): _ZETA3, s.zeta(5): _ZETA5})
        polynomial = s.Poly(s.expand(formal), _PHASE, _ZETA3, _ZETA5)
        for monomial, value in polynomial.terms():
            if value.is_Rational is not True:
                raise UnsupportedExpressionError("Named-function constants exceed the exact basis")
            if value:
                features[word, monomial] = value
    return features


@lru_cache(maxsize=5)
def _basis(weight: int):
    candidates = tuple(
        sorted(
            set(_monomials(weight)),
            key=lambda value: (s.count_ops(value), s.default_sort_key(value)),
        )
    )
    columns = [_features(candidate) for candidate in candidates]
    rows = tuple(sorted(set().union(*columns)))
    matrix = s.ImmutableMatrix([[column.get(row, 0) for column in columns] for row in rows])
    return candidates, rows, matrix


@lru_cache(maxsize=180)
def _normal_form(n: int, p: int, mode: int, boundary: Boundary) -> s.Expr:
    argument = _ARGUMENTS[mode]
    source = s.polylog(n + 1, argument) if p == 1 else NielsenPolylog(n, p, argument)
    target = _features(source, boundary)
    candidates, rows, matrix = _basis(n + p)
    if set(target).difference(rows):
        raise UnsupportedExpressionError("Polylogarithm exceeds the named-function basis")
    try:
        solution = matrix.gauss_jordan_solve(
            s.ImmutableMatrix([target.get(row, 0) for row in rows])
        )
    except ValueError as exc:
        raise UnsupportedExpressionError(
            "Polylogarithm has no proved named-function reduction"
        ) from exc
    coefficients, parameters = solution[0], solution[1]
    if parameters.rows:
        raise UnsupportedExpressionError("Named-function basis is unexpectedly dependent")
    return s.Add(
        *(coefficients[index, 0] * candidate for index, candidate in enumerate(candidates))
    )


def reduce_polylogs(
    expression: s.Expr,
    z: s.Symbol,
    *,
    domain: Interval,
    boundary: Boundary | None = None,
) -> s.Expr:
    """Reduce classical and Nielsen polylogarithms through weight five.

    On domain='0<z<1', reduce atoms at z, 1-z, z/(z-1), 1-1/z, 1/z and
    1/(1-z) to logs, classical polylogarithms and the single Nielsen function
    S_{3,2}(z). In particular S_{2,2} becomes classical and S_{2,3} reduces to
    that weight-five basis. The coefficients are generated by exact shuffle
    identities and rational linear algebra, with no stored reduction table.

    For arguments above one, boundary='upper' or 'lower' specifies the side
    of each original polylogarithm argument. Reduced atoms have real arguments
    below one; explicit i*pi terms carry the cut data.
    Higher weights and unrecognized arguments remain unchanged. Nielsen
    constants at zero and +/-1 use reduce_nielsen's exact endpoint reductions.

    Args:
        expression: Exact expression containing classical or Nielsen polylogs.
        z: Common real interior variable, a commutative SymPy Symbol.
        domain: Explicit declaration '0<z<1'.
        boundary: Required side for transformed arguments above one.
    """
    if not isinstance(expression, s.Expr) or not isinstance(z, s.Symbol):
        raise TypeError("Use an exact SymPy expression and symbolic target argument")
    if z.is_commutative is not True:
        raise TypeError("The target argument must be commutative")
    if domain != "0<z<1":
        raise UnsupportedExpressionError("Explicit domain='0<z<1' is required")
    if boundary not in (None, "upper", "lower"):
        raise ValueError("boundary must be upper, lower, or omitted")
    if expression.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError("Polylogarithm reduction requires exact finite input")
    if s.count_ops(expression) > 4096:
        raise UnsupportedExpressionError("Polylogarithm expression exceeds the operation budget")
    expression = reduce_nielsen(expression)
    replacements = {}
    arguments = tuple(cast(s.Expr, argument.xreplace({_Z: z})) for argument in _ARGUMENTS)
    for atom in expression.atoms(s.polylog, NielsenPolylog):
        argument = cast(s.Expr, atom.args[-1])
        if not argument.has(z):
            continue
        if isinstance(atom, NielsenPolylog):
            n, p = int(atom.n), int(atom.p)
        else:
            order = atom.args[0]
            if not isinstance(order, s.Integer) or not 1 <= order <= 5:
                continue
            n, p = int(order) - 1, 1
        if n + p > 5:
            continue
        mode = next(
            (index for index, target in enumerate(arguments) if s.cancel(argument - target) == 0),
            None,
        )
        if mode is None:
            continue
        if mode in (4, 5) and boundary is None:
            raise UnsupportedExpressionError("Arguments above one require an explicit boundary")
        replacements[atom] = _normal_form(n, p, mode, boundary or "upper").xreplace({_Z: z})
    return cast(s.Expr, expression.xreplace(replacements))
