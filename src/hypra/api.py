"""Public expansion API with exact input and checked Laurent output."""

from dataclasses import asdict, dataclass
from typing import Any

import sympy as s

from .errors import IncompleteExpansionError, UnsupportedExpressionError
from .planning import laurent_coefficient, plan_expansion


@dataclass(frozen=True)
class ExpansionResult:
    """Exact Laurent polynomial plus its remainder and reproducibility record.

    ``expression`` is a SymPy expression including O(epsilon**(through+1)).
    ``polynomial`` omits the remainder, for coefficient algebra. The result's
    ``_sympy_`` hook also allows ``sympy.sympify(result)``.
    """

    polynomial: s.Expr
    epsilon: s.Symbol
    through: int
    provenance: dict[str, Any]

    @property
    def expression(self):
        return self.polynomial + s.Order(self.epsilon ** (self.through + 1), (self.epsilon, 0))

    def coefficient(self, power: int):
        if not isinstance(power, int) or isinstance(power, bool):
            raise TypeError("power must be an integer")
        if power > self.through:
            raise ValueError("Coefficient is beyond the certified truncation order")
        return laurent_coefficient(self.polynomial, self.epsilon, power)

    def latex(self, *, include_order: bool = True) -> str:
        """Render exact mathematics, including symbolic HPL words, as LaTeX."""
        return s.latex(self.expression if include_order else self.polynomial)

    def _sympy_(self):
        return self.expression

    def to_dict(self):
        from .serialization import encode_expr

        return {
            "polynomial": encode_expr(self.polynomial),
            "epsilon": encode_expr(self.epsilon),
            "through": self.through,
            "remainder_order": self.through + 1,
            "latex": self.latex(),
            "provenance": self.provenance,
        }


def backend_info(config=None, *, backend="native"):
    """Report native backend capabilities, or explicitly inspect external Wolfram."""
    if backend == "native":
        from .native import NativeBackend

        active = NativeBackend(config)
    elif backend == "wolfram":
        from .backend import WolframBackend

        active = WolframBackend(config)
    else:
        raise ValueError("backend must be 'native' or 'wolfram'")
    info = active.info()
    return info.as_dict() if hasattr(info, "as_dict") else asdict(info)


def _checked_polynomial(expr, epsilon, through, lower_bound):
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


def expand_epsilon(expr, epsilon, through=3, assumptions=None, *, config=None, backend=None):
    """Expand the complete expression through the inclusive Laurent order.

    The default backend is native Python/SymPy and supports documented integer and half-integer
    parameter families. Every prefactor is included; atoms are expanded deeper
    when poles require it. Unsupported cases fail closed without fallback.
    Set backend="wolfram" explicitly to use the experimental external bridge.

    Parameters
    ----------
    expr : sympy.Expr
        Exact symbolic input (no strings or binary floating-point atoms).
    epsilon : sympy.Symbol
        Regulator expanded at zero.
    through : int
        Inclusive order, e.g. 3 includes all negative powers and epsilon**3.
    assumptions : SymPy Boolean, optional
        Explicit relational assumptions, e.g. And(w > 0, w < 1).
        Symbol assumptions on inputs are also preserved. This API does not
        exchange integrals/limits or change branches by power simplification.
    config : NativeConfig or BackendConfig, optional
        Configuration matching the selected backend.
    backend : "native", "wolfram" or backend object, optional
        Native is the default. Custom objects are marked as injected in
        provenance; there is no silent substitution between implementations.
    """
    from .backend import WolframBackend
    from .serialization import decode_expr, encode_expr

    if not isinstance(expr, s.Expr):
        raise TypeError("expr must be an exact SymPy expression, not a string")
    if not isinstance(epsilon, s.Symbol):
        raise TypeError("epsilon must be a SymPy Symbol")
    if assumptions is None:
        assumptions = s.true
    if not isinstance(assumptions, s.logic.boolalg.Boolean):
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
    declared = {}
    for symbol in expr.free_symbols | assumptions.free_symbols | {epsilon}:
        if symbol.name in declared and declared[symbol.name] != symbol:
            raise UnsupportedExpressionError(
                "Symbols with the same name and different assumptions are ambiguous"
            )
        declared[symbol.name] = symbol
    original_ast = expression_ast
    if backend is None or backend == "native":
        from .endpoints import rewrite_endpoints

        expr = rewrite_endpoints(expr, epsilon, assumptions)
        expression_ast = encode_expr(expr)
    plan = plan_expansion(expr, epsilon, through, assumptions)
    injected = backend is not None and not isinstance(backend, str)
    if backend is None or backend == "native":
        from .native import NativeBackend

        active = NativeBackend(config)
    elif backend == "wolfram":
        active = WolframBackend(config)
    elif isinstance(backend, str):
        raise ValueError("Unknown backend; choose 'native' or 'wolfram'")
    else:
        active = backend
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
    if any(symbol.name not in declared for symbol in polynomial.free_symbols):
        raise IncompleteExpansionError("Backend introduced undeclared symbols")
    polynomial = polynomial.xreplace(
        {symbol: declared[symbol.name] for symbol in polynomial.free_symbols}
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
        "branch_convention": "Principal branches; no automatic analytic continuation or integral/limit interchange",
        "certificate_scope": "Conditional on backend correctness and the stated input domain",
    }
    return ExpansionResult(polynomial, epsilon, through, provenance)
