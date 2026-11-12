"""Conservative Laurent-order planning; no numerical inference of pole orders."""

from dataclasses import dataclass

import sympy as s

from .errors import UnsupportedExpressionError


@dataclass(frozen=True)
class ExpansionPlan:
    through: int
    pole_bound: int
    internal_order: int
    valuation_lower_bound: int


def _facts(assumptions):
    clauses = assumptions.args if isinstance(assumptions, s.And) else (assumptions,)
    facts = []
    for clause in clauses:
        if isinstance(clause, s.StrictGreaterThan):
            facts.append(s.Q.positive(clause.lhs - clause.rhs))
        elif isinstance(clause, s.StrictLessThan):
            facts.append(s.Q.positive(clause.rhs - clause.lhs))
        elif isinstance(clause, s.Unequality):
            facts.append(s.Q.nonzero(clause.lhs - clause.rhs))
        else:
            facts.append(clause)
    return s.And(*facts)


def _positive(expr, assumptions):
    return s.ask(s.Q.positive(expr), _facts(assumptions)) is True


def _affine_rational(expr, epsilon):
    try:
        poly = s.Poly(expr, epsilon)
    except s.PolynomialError as exc:
        raise UnsupportedExpressionError(
            "Hypergeometric parameters must be rational affine in epsilon"
        ) from exc
    if poly.degree() > 1 or not all(c.is_Rational for c in poly.all_coeffs()):
        raise UnsupportedExpressionError(
            "Hypergeometric parameters must be rational affine in epsilon"
        )
    base = poly.nth(0)
    if base.q not in (1, 2):
        raise UnsupportedExpressionError("Only integer/half-integer parameter bases are supported")
    return base, poly.nth(1)


def _exact_rational_valuation(expr, epsilon, assumptions):
    """Order of a rational polynomial quotient when leading coefficients are nonzero."""
    try:
        num, den = s.cancel(expr).as_numer_denom()
        pn, pd = s.Poly(num, epsilon), s.Poly(den, epsilon)
    except s.PolynomialError as exc:
        raise UnsupportedExpressionError("Cannot establish denominator order exactly") from exc
    if pn.is_zero:
        raise UnsupportedExpressionError(
            "A zero denominator or inverse-zero expression is undefined"
        )
    orders = []
    for p in (pn, pd):
        k = min(m[0] for m, c in p.terms() if c != 0)
        c = p.nth(k)
        if c.is_nonzero is not True and s.ask(s.Q.nonzero(c), _facts(assumptions)) is not True:
            raise UnsupportedExpressionError(
                f"Nonzero leading coefficient {c} requires an explicit nonzero assumption"
            )
        orders.append(k)
    return orders[0] - orders[1]


def plan_expansion(expr, epsilon, through, assumptions=s.true):
    """Bound every omitted hypergeometric term after multiplying all prefactors.

    If the total negative valuation bound is p, request each hypergeometric
    atom through N+p. This may over-expand but cannot lose coefficients through
    N under the supported affine-parameter/meromorphic-input restrictions.
    Unknown orders and essential singularities are rejected, not guessed.
    """
    if not isinstance(through, int) or isinstance(through, bool):
        raise TypeError("through must be an integer inclusive Laurent order")
    negative_budget = 0
    # Validate even regulator-independent subexpressions, which valuation
    # recursion legitimately treats as constants but must not define at poles.
    if any(
        not b.has(epsilon) and b.is_integer is True and b.is_nonpositive is True
        for atom in expr.atoms(s.hyper) for b in atom.bq
    ):
        raise UnsupportedExpressionError("Fixed singular lower hypergeometric parameter")

    def visit(node):
        nonlocal negative_budget
        if not node.has(epsilon):
            return 0
        if node == epsilon:
            return 1
        if isinstance(node, s.Add):
            return min(visit(a) for a in node.args)
        if isinstance(node, s.Mul):
            return sum(visit(a) for a in node.args)
        if isinstance(node, s.Pow):
            base, power = node.args
            if power.is_Integer:
                if power >= 0:
                    # Each factor in a repeated product may multiply an omitted term.
                    before = negative_budget
                    value = visit(base)
                    negative_budget = before + (negative_budget - before) * int(power)
                    return int(power) * value
                if base.func in (s.gamma, s.exp):
                    before = negative_budget
                    val = visit(base) * int(power)
                    negative_budget = before
                else:
                    val = _exact_rational_valuation(base, epsilon, assumptions) * int(power)
                negative_budget += max(0, -val)
                return val
            if power.has(epsilon):
                if base.has(epsilon):
                    raise UnsupportedExpressionError(
                        "Regulator-dependent base and exponent are unsupported"
                    )
                # exp(power*log(base)) is analytic in epsilon for nonzero fixed base,
                # but keeping a real positive base avoids unqualified branch changes.
                if not _positive(base, assumptions):
                    raise UnsupportedExpressionError(
                        f"Regulator-dependent powers require a positive base assumption: {base}"
                    )
                if visit(power) < 0:
                    raise UnsupportedExpressionError("Essential singularities are unsupported")
                return 0
            raise UnsupportedExpressionError(
                "Fractional epsilon-dependent powers are not Laurent series"
            )
        if node.func == s.gamma:
            arg = node.args[0]
            try:
                p = s.Poly(arg, epsilon)
            except s.PolynomialError as exc:
                raise UnsupportedExpressionError(
                    "Gamma argument must be polynomial in epsilon"
                ) from exc
            c = p.nth(0)
            if c.is_Rational:
                if c.is_integer and c <= 0:
                    powers = [k[0] for k, value in p.terms() if k[0] > 0 and value != 0]
                    if not powers:
                        raise UnsupportedExpressionError("Gamma has a fixed pole")
                    m = min(powers)
                    coeff = p.nth(m)
                    if (
                        coeff.is_nonzero is not True
                        and s.ask(s.Q.nonzero(coeff), _facts(assumptions)) is not True
                    ):
                        raise UnsupportedExpressionError(
                            "Gamma pole slope must be provably nonzero"
                        )
                    negative_budget += m
                    return -m
                return 0
            if not _positive(c, assumptions):
                raise UnsupportedExpressionError(
                    "Symbolic Gamma base requires a positive assumption"
                )
            return 0
        if node.func == s.hyper:
            if len(node.ap) != len(node.bq) + 1:
                raise UnsupportedExpressionError("Initial backend scope is pF(p-1)")
            if node.argument.has(epsilon):
                raise UnsupportedExpressionError(
                    "Hypergeometric argument must be epsilon-independent"
                )
            for a in node.ap:
                _affine_rational(a, epsilon)
            poles = 0
            for b in node.bq:
                c, slope = _affine_rational(b, epsilon)
                if c.is_integer and c <= 0:
                    if slope == 0:
                        raise UnsupportedExpressionError(
                            "Fixed singular lower hypergeometric parameter"
                        )
                    poles += 1
            negative_budget += poles
            return -poles
        if node.func == s.exp:
            if visit(node.args[0]) < 0:
                raise UnsupportedExpressionError("Essential singularities are unsupported")
            return 0
        if node.func == s.log:
            arg = node.args[0]
            if not _positive(arg.subs(epsilon, 0), assumptions):
                raise UnsupportedExpressionError(
                    "Log argument at epsilon=0 must be provably positive"
                )
            if visit(arg) < 0:
                raise UnsupportedExpressionError(
                    "Logarithmic epsilon singularities are not Laurent series"
                )
            return 0
        raise UnsupportedExpressionError(f"Cannot establish Laurent order for {node.func}")

    lower = visit(expr)
    # At least Taylor order zero is needed even for a negative requested order.
    internal = max(0, through + negative_budget)
    return ExpansionPlan(through, negative_budget, internal, lower)


def laurent_coefficient(expr, epsilon, power):
    """Extract a Laurent coefficient without assuming epsilon is factored visually.

    SymPy expand may put epsilon inside a rational denominator; its Expr.coeff
    then misses poles. as_coeff_exponent separates those powers exactly.
    """
    result = s.S.Zero
    for term in s.Add.make_args(s.expand(expr)):
        if term == 0:
            continue
        coefficient, exponent = term.as_coeff_exponent(epsilon)
        if coefficient.has(epsilon) or not exponent.is_Integer:
            raise UnsupportedExpressionError("Expression is not an explicit Laurent polynomial")
        if exponent == power:
            result += coefficient
    return result
