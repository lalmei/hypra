"""Public expansion API with exact input and checked Laurent output."""

from dataclasses import asdict, dataclass
from typing import Any, Literal, Protocol, cast

import sympy as s
from sympy.logic.boolalg import Boolean

from .config import BackendConfig
from .errors import IncompleteExpansionError, UnsupportedExpressionError
from .native import NativeBackend, NativeConfig
from .native.planning import laurent_coefficient, plan_expansion

type BackendName = Literal["native", "wolfram"]


class ExpansionBackend(Protocol):
    """Structural interface for an injected expansion backend."""

    def expand(self, request: dict[str, Any], /) -> dict[str, Any]:
        """Return a certified expansion response for the exact AST request."""
        ...


@dataclass(frozen=True)
class ExpansionResult:
    """Exact Laurent polynomial plus its remainder and reproducibility record.

    ``sympy.sympify(result)`` includes the order remainder.

    Attributes:
        polynomial: Exact Laurent polynomial used for coefficient algebra.
        epsilon: Regulator expanded about zero.
        through: Highest certified Laurent power.
        provenance: Input ASTs, plan, backend details, and branch conventions.
    """

    polynomial: s.Expr
    epsilon: s.Symbol
    through: int
    provenance: dict[str, Any]

    @property
    def expression(self) -> s.Expr:
        """Exact Laurent expansion including the SymPy order remainder."""
        return self.polynomial + s.Order(self.epsilon ** (self.through + 1), (self.epsilon, 0))

    def coefficient(self, power: int) -> s.Expr:
        """Return the coefficient of the requested Laurent power.

        Args:
            power: Integer exponent of the regulator.

        Raises:
            TypeError: The exponent is not an integer.
            ValueError: The exponent exceeds the certified order.
        """
        if not isinstance(power, int) or isinstance(power, bool):
            raise TypeError("power must be an integer")
        if power > self.through:
            raise ValueError("Coefficient is beyond the certified truncation order")
        return laurent_coefficient(self.polynomial, self.epsilon, power)

    def latex(self, *, include_order: bool = True) -> str:
        """Render exact mathematics, including symbolic HPL words, as LaTeX.

        Args:
            include_order: Include the SymPy order remainder.
        """
        return s.latex(self.expression if include_order else self.polynomial)

    def _sympy_(self) -> s.Expr:
        return self.expression

    def to_dict(self) -> dict[str, Any]:
        """Return the exact JSON-compatible ASTs, LaTeX, and provenance."""
        from .serialization import encode_expr

        return {
            "polynomial": encode_expr(self.polynomial),
            "epsilon": encode_expr(self.epsilon),
            "through": self.through,
            "remainder_order": self.through + 1,
            "latex": self.latex(),
            "provenance": self.provenance,
        }


def backend_info(
    config: NativeConfig | BackendConfig | None = None,
    *,
    backend: BackendName = "native",
) -> dict[str, Any]:
    """Report native backend capabilities or inspect an external Wolfram installation.

    Args:
        config: Configuration for the selected backend.
        backend: Engine to inspect.

    Returns:
        A JSON-compatible capability record including engine availability.

    Raises:
        ValueError: The backend name is unknown.
    """
    if backend == "native":
        info = NativeBackend(cast(NativeConfig | None, config)).info()
    elif backend == "wolfram":
        from .backend import WolframBackend

        info = WolframBackend(cast(BackendConfig | None, config)).info()
    else:
        raise ValueError("backend must be 'native' or 'wolfram'")
    return info.as_dict() if hasattr(info, "as_dict") else asdict(info)


def _checked_polynomial(
    expr: s.Expr,
    epsilon: s.Symbol,
    through: int,
    lower_bound: int,
) -> s.Expr:
    if expr.has(
        s.Float, s.nan, s.zoo, s.oo, -s.oo, s.Derivative, s.Subs, s.Integral, s.Sum, s.hyper
    ):
        raise IncompleteExpansionError(
            "Backend returned non-exact, unresolved, or non-finite output"
        )
    expanded = s.expand(expr)
    for term in s.Add.make_args(expanded):
        if term == 0:
            continue
        coefficient, power = term.as_coeff_exponent(epsilon)
        if coefficient.has(epsilon) or not power.is_Integer:
            raise IncompleteExpansionError("Backend output is not a Laurent polynomial in epsilon")
        if power > through or (term != 0 and power < lower_bound):
            raise IncompleteExpansionError(
                "Backend output contradicts the requested Laurent bounds"
            )
    return expanded


def expand_epsilon(
    expr: s.Expr,
    epsilon: s.Symbol,
    through: int = 3,
    assumptions: Boolean | None = None,
    *,
    config: NativeConfig | BackendConfig | None = None,
    backend: BackendName | ExpansionBackend | None = None,
) -> ExpansionResult:
    """Expand the complete expression through the inclusive Laurent order.

    The default backend is native Python/SymPy and supports documented integer and half-integer
    parameter families. Every prefactor is included; atoms are expanded deeper
    when poles require it. Unsupported cases fail closed without fallback.
    Set backend="wolfram" explicitly to use the experimental external bridge.

    Args:
        expr: Exact symbolic input without floating-point atoms.
        epsilon: Regulator expanded at zero.
        through: Inclusive order. For example, 3 includes all poles and epsilon**3.
        assumptions: Explicit relational assumptions, such as `And(w > 0, w < 1)`.
            Assumptions attached to input symbols are preserved.
        config: Settings for the selected backend. Native settings can opt in
            to meromorphic parameter continuation at z=1.
        backend: Native, Wolfram, or a custom engine implementing `expand`.
            The result records a custom engine as injected.

    Returns:
        Exact Laurent polynomial, order remainder, and provenance.

    Raises:
        TypeError: The input, regulator, assumptions, or configuration has the wrong type.
        ValueError: The backend name is unknown.
        UnsupportedExpressionError: The expression or its domain is unsupported.
        IncompleteExpansionError: The backend could not certify a complete exact result.
    """
    from .backend import WolframBackend
    from .serialization import decode_expr, encode_expr

    if not isinstance(expr, s.Expr):
        raise TypeError("expr must be an exact SymPy expression, not a string")
    if not isinstance(epsilon, s.Symbol):
        raise TypeError("epsilon must be a SymPy Symbol")
    if assumptions is None:
        assumptions = s.true
    if not isinstance(assumptions, Boolean):
        raise TypeError("assumptions must be a SymPy Boolean expression")
    if s.simplify(assumptions) is s.false:
        raise UnsupportedExpressionError("Contradictory assumptions are not a valid domain")
    try:
        if s.satisfiable(assumptions, use_lra_theory=True) is False:
            raise UnsupportedExpressionError("Contradictory assumptions are not a valid domain")
    except (AssertionError, ValueError, NotImplementedError) as exc:
        if isinstance(exc, UnsupportedExpressionError):
            raise
        # Isolated complex-symbol nonzero guards need no ordering theory.
        # Remove only Ne(symbol,0) guards whose symbols occur nowhere else in
        # the assumption formula; all remaining domain relations still pass LRA.
        clauses = list(assumptions.args) if isinstance(assumptions, s.And) else [assumptions]
        independent = []
        for clause in clauses:
            if not isinstance(clause, s.Unequality):
                continue
            symbol = clause.lhs if clause.rhs == 0 else clause.rhs if clause.lhs == 0 else None
            if isinstance(symbol, s.Symbol) and all(
                other == clause or not other.has(symbol) for other in clauses
            ):
                independent.append(clause)
        if not independent:
            raise UnsupportedExpressionError(
                "Assumption consistency is unsupported; use linear real-domain relations"
            ) from exc
        remaining = s.And(*(clause for clause in clauses if clause not in independent))
        try:
            if s.satisfiable(remaining, use_lra_theory=True) is False:
                raise UnsupportedExpressionError("Contradictory assumptions are not a valid domain")
        except (AssertionError, ValueError, NotImplementedError) as remainder_error:
            if isinstance(remainder_error, UnsupportedExpressionError):
                raise
            raise UnsupportedExpressionError(
                "Assumption consistency is unsupported beyond isolated nonzero guards"
            ) from remainder_error
    if expr.has(s.Float, s.nan, s.zoo, s.oo, -s.oo):
        raise UnsupportedExpressionError(
            "Input must be exact and finite; use Rational instead of Float"
        )
    # Encode before invoking anything external: this checks all input heads,
    # identifiers, exact constants and assumption syntax against a whitelist.
    expression_ast, epsilon_ast = encode_expr(expr), encode_expr(epsilon)
    assumption_ast = encode_expr(assumptions)
    declared: dict[str, s.Symbol] = {}
    input_symbols = cast(set[s.Symbol], expr.free_symbols | assumptions.free_symbols | {epsilon})
    for symbol in input_symbols:
        if symbol.name in declared and declared[symbol.name] != symbol:
            raise UnsupportedExpressionError(
                "Symbols with the same name and different assumptions are ambiguous"
            )
        declared[symbol.name] = symbol
    original_ast = expression_ast
    native: NativeBackend | None = None
    if backend is None or backend == "native":
        from .native.endpoints import rewrite_endpoints

        native = NativeBackend(cast(NativeConfig | None, config))
        expr = rewrite_endpoints(
            expr,
            epsilon,
            assumptions,
            meromorphic=native.config.unit_endpoint == "meromorphic",
            excess_shift=native.config.unit_excess_shift,
        )
        expression_ast = encode_expr(expr)
    plan = plan_expansion(expr, epsilon, through, assumptions)
    active: ExpansionBackend
    if backend == "wolfram":
        active = WolframBackend(cast(BackendConfig | None, config))
        injected = False
    elif backend is not None and not isinstance(backend, str):
        active = backend
        injected = True
    elif native is not None:
        active = native
        injected = False
    else:
        raise ValueError("Unknown backend; choose 'native' or 'wolfram'")
    request = {
        "protocol": 1,
        "operation": "expand",
        "expression": expression_ast,
        "epsilon": epsilon_ast,
        "through": through,
        "internal_order": plan.internal_order,
        "assumptions": assumption_ast,
    }
    response = active.expand(request)
    if not isinstance(response, dict):
        raise IncompleteExpansionError("Backend result must be a structured object")
    if response.get("complete") is not True:
        raise IncompleteExpansionError("Backend marked the expansion incomplete")
    if (
        type(response.get("achieved_through")) is not int
        or response["achieved_through"] != through
        or type(response.get("remainder_order")) is not int
        or response["remainder_order"] != through + 1
    ):
        raise IncompleteExpansionError("Backend did not certify the requested inclusive order")
    try:
        polynomial = decode_expr(response["result"])
    except (KeyError, ValueError, TypeError) as exc:
        raise IncompleteExpansionError("Backend returned an invalid symbolic result") from exc
    if not isinstance(polynomial, s.Expr):
        raise IncompleteExpansionError("Backend result must be a scalar SymPy expression")
    output_symbols = cast(set[s.Symbol], polynomial.free_symbols)
    if any(symbol.name not in declared for symbol in output_symbols):
        raise IncompleteExpansionError("Backend introduced undeclared symbols")
    polynomial = cast(
        s.Expr,
        polynomial.xreplace({symbol: declared[symbol.name] for symbol in output_symbols}),
    )
    polynomial = _checked_polynomial(polynomial, epsilon, through, plan.valuation_lower_bound)
    backend_record = response.get("backend", {})
    if not isinstance(backend_record, dict):
        raise IncompleteExpansionError("Backend provenance must be a structured object")
    provenance = {
        "hypra_version": "0.1.0a0",
        "sympy_version": s.__version__,
        "input": original_ast,
        "planned_input": expression_ast,
        "epsilon": epsilon_ast,
        "assumptions": assumption_ast,
        "plan": asdict(plan),
        "backend": backend_record,
        "backend_fingerprint": response.get("provenance", {}),
        "cache_hit": response.get("cache_hit", False),
        "order_basis": response.get("order_basis"),
        "injected_backend": injected,
        "branch_convention": "Principal branches",
        "unit_endpoint": native.config.unit_endpoint if native is not None else None,
        "unit_excess_shift": native.config.unit_excess_shift if native is not None else None,
        "certificate_scope": "Conditional on backend correctness and the stated input domain",
    }
    return ExpansionResult(polynomial, epsilon, through, provenance)
