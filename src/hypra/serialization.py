"""Strict data-only JSON AST transport between SymPy and Wolfram.

Never uses eval, sympify on strings, parse_expr, or executable source text.
Numeric atoms are decimal strings to avoid JSON integer-range/float loss.
Truncation orders and analytic assumptions belong to the request/result
protocol; ``Order`` and unevaluated integrals/derivatives are not accepted.
"""

from __future__ import annotations

import re
from collections.abc import Mapping

import sympy as sp
from sympy.core.assumptions import _assume_defined
from sympy.core.relational import Relational
from sympy.logic.boolalg import Boolean

from .symbolic import HPL, MZV, EulerSum


class SerializationError(ValueError):
    """An expression or AST cannot be transported exactly and safely."""


MAX_DEPTH = 128
MAX_NODES = 100000
MAX_DIGITS = 10000
_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,255}\Z")
_INTEGER = re.compile(r"(?:0|-[1-9][0-9]*|[1-9][0-9]*)\Z")
_CONSTANTS = {
    sp.pi: "Pi",
    sp.E: "E",
    sp.I: "ImaginaryUnit",
    sp.EulerGamma: "EulerGamma",
    sp.Catalan: "Catalan",
    sp.true: "True",
    sp.false: "False",
}
_CONSTANTS_BY_HEAD = {head: value for value, head in _CONSTANTS.items()}
_FUNCTIONS = {
    sp.Add: "Add",
    sp.Mul: "Mul",
    sp.Pow: "Pow",
    sp.gamma: "Gamma",
    sp.log: "Log",
    sp.polylog: "PolyLog",
    sp.zeta: "Zeta",
    sp.exp: "Exp",
    sp.polygamma: "PolyGamma",
    sp.And: "And",
    sp.Or: "Or",
    sp.Not: "Not",
    sp.Eq: "Equal",
    sp.Ne: "Unequal",
    sp.Lt: "Less",
    sp.Le: "LessEqual",
    sp.Gt: "Greater",
    sp.Ge: "GreaterEqual",
}
_FUNCTIONS_BY_HEAD = {head: fn for fn, head in _FUNCTIONS.items()}
_ARITY = {
    "Pow": (2,),
    "Gamma": (1,),
    "Log": (1, 2),
    "PolyLog": (2,),
    "Zeta": (1, 2),
    "Exp": (1,),
    "PolyGamma": (2,),
    "Not": (1,),
    "Equal": (2,),
    "Unequal": (2,),
    "Less": (2,),
    "LessEqual": (2,),
    "Greater": (2,),
    "GreaterEqual": (2,),
    "HypergeometricPFQ": (3,),
    "HPL": (2,),
    "MZV": (1,),
    "EulerSum": (2,),
}
_RELATIONS = {"Equal", "Unequal", "Less", "LessEqual", "Greater", "GreaterEqual"}
_LOGIC = {"And", "Or", "Not"}


def _name(name):
    if not isinstance(name, str) or not _NAME.fullmatch(name):
        raise SerializationError("Symbol names must be ASCII letters followed by letters/digits/_")
    return name


def _assumptions(facts):
    if not isinstance(facts, dict):
        raise SerializationError("Symbol assumptions must be a dictionary")
    if any(key not in _assume_defined or type(value) is not bool for key, value in facts.items()):
        raise SerializationError("Symbol assumptions must contain known Boolean SymPy facts")
    return facts


def _number(value):
    if not isinstance(value, str) or len(value) > MAX_DIGITS or not _INTEGER.fullmatch(value):
        raise SerializationError("Exact integer must be a canonical decimal string")
    try:
        return int(value)
    except ValueError as exc:
        raise SerializationError(
            "Exact integer exceeds this Python runtime's supported size"
        ) from exc


def encode_expr(expression):
    """Encode a supported exact SymPy expression or Boolean condition.

    Python integers and booleans are accepted. Strings, floats, integrals,
    derivatives, unevaluated orders and all unknown function heads are rejected.
    Symbol names deliberately exclude Wolfram contexts and source syntax.
    """
    budget = [MAX_NODES]

    def encode(expr, depth):
        budget[0] -= 1
        if depth > MAX_DEPTH or budget[0] < 0:
            raise SerializationError("Expression exceeds AST size/depth limit")
        if type(expr) is bool:
            expr = sp.true if expr else sp.false
        elif type(expr) is int:
            expr = sp.Integer(expr)
        if not isinstance(expr, sp.Basic):
            raise SerializationError(f"Unsupported input type: {type(expr).__name__}")
        if isinstance(expr, sp.Float):
            raise SerializationError("Floating-point input is forbidden; use exact Rational values")
        if isinstance(expr, sp.Integer):
            return {"head": "Integer", "value": str(expr)}
        if isinstance(expr, sp.Rational):
            return {"head": "Rational", "numerator": str(expr.p), "denominator": str(expr.q)}
        if isinstance(expr, sp.Symbol):
            if type(expr) is not sp.Symbol or expr.is_commutative is False:
                raise SerializationError("Only ordinary commutative Symbols are supported")
            # Use original (not inferred) assumptions so reconstruction is reproducible.
            facts = {k: v for k, v in expr._assumptions_orig.items() if v is not None}
            return {"head": "Symbol", "name": _name(expr.name), "assumptions": _assumptions(facts)}
        if expr in _CONSTANTS:
            return {"head": _CONSTANTS[expr]}
        if isinstance(expr, sp.Tuple):
            head = "List"
        elif expr.func == sp.hyper:
            head = "HypergeometricPFQ"
        elif expr.func == HPL:
            head = "HPL"
        elif expr.func == EulerSum:
            head = "EulerSum"
        elif expr.func == MZV:
            head = "MZV"
        else:
            head = _FUNCTIONS.get(expr.func)
        if head is None:
            raise SerializationError(f"Unsupported symbolic head: {expr.func.__name__}")
        return {"head": head, "args": [encode(arg, depth + 1) for arg in expr.args]}

    return encode(expression, 0)


def decode_expr(payload):
    """Reconstruct a whitelisted data-only AST, rejecting malformed fields.

    This API accepts already parsed JSON values. JSON parsers should also reject
    duplicate object keys before calling it. No node can contain executable code.
    """
    budget = [MAX_NODES]

    def decode(node, depth):
        budget[0] -= 1
        if depth > MAX_DEPTH or budget[0] < 0:
            raise SerializationError("AST exceeds size/depth limit")
        if not isinstance(node, Mapping) or not isinstance(node.get("head"), str):
            raise SerializationError("Every AST node must be an object with a string head")
        head = node["head"]
        if head == "Integer":
            _keys(node, {"head", "value"})
            return sp.Integer(_number(node["value"]))
        if head == "Rational":
            _keys(node, {"head", "numerator", "denominator"})
            numerator, denominator = _number(node["numerator"]), _number(node["denominator"])
            if denominator <= 0:
                raise SerializationError("Rational denominator must be strictly positive")
            return sp.Rational(numerator, denominator)
        if head == "Symbol":
            _keys(node, {"head", "name"}, {"assumptions"})
            try:
                facts = _assumptions(node.get("assumptions", {}))
                if facts.get("commutative") is False:
                    raise SerializationError("Noncommutative Symbols are not supported")
                return sp.Symbol(_name(node["name"]), **facts)
            except (TypeError, ValueError) as exc:
                raise SerializationError("Invalid or inconsistent symbol assumptions") from exc
        if head in _CONSTANTS_BY_HEAD:
            _keys(node, {"head"})
            return _CONSTANTS_BY_HEAD[head]
        if head not in _FUNCTIONS_BY_HEAD and head not in {
            "List",
            "HypergeometricPFQ",
            "HPL",
            "MZV",
            "EulerSum",
        }:
            raise SerializationError(f"Unsupported AST head: {head}")
        _keys(node, {"head", "args"})
        if not isinstance(node["args"], list):
            raise SerializationError("AST args must be a JSON array")
        if head in _ARITY and len(node["args"]) not in _ARITY[head]:
            raise SerializationError(f"Wrong argument count for {head}")
        args = [decode(arg, depth + 1) for arg in node["args"]]
        if head == "List":
            return sp.Tuple(*args)
        if head == "HypergeometricPFQ":
            if not isinstance(args[0], sp.Tuple) or not isinstance(args[1], sp.Tuple):
                raise SerializationError("Hypergeometric upper/lower parameters must be Lists")
            _scalar(args[2])
            for arg in (*args[0], *args[1]):
                _scalar(arg)
            return sp.hyper(args[0], args[1], args[2], evaluate=False)
        if head == "HPL":
            if not isinstance(args[0], sp.Tuple):
                raise SerializationError("HPL word must be a List")
            _scalar(args[1])
            try:
                return HPL(args[0], args[1])
            except (TypeError, ValueError) as exc:
                raise SerializationError(
                    "HPL word must contain expanded letters -1, 0, or 1"
                ) from exc
        if head == "EulerSum":
            try:
                return EulerSum(*args)
            except (TypeError, ValueError) as exc:
                raise SerializationError("Invalid or divergent EulerSum") from exc
        if head == "MZV":
            if not isinstance(args[0], sp.Tuple):
                raise SerializationError("MZV indices must be a List")
            try:
                return MZV(args[0])
            except (TypeError, ValueError) as exc:
                raise SerializationError(
                    "MZV requires positive indices with first index >= 2"
                ) from exc
        if head in _LOGIC:
            if any(not isinstance(arg, Boolean) for arg in args):
                raise SerializationError("Logical operators require Boolean arguments")
        else:
            for arg in args:
                _scalar(arg)
        try:
            # Unevaluated relations avoid dropping assumptions during transport.
            result = (
                _FUNCTIONS_BY_HEAD[head](*args, evaluate=False)
                if head in _RELATIONS
                else (_FUNCTIONS_BY_HEAD[head](*args))
            )
        except (TypeError, ValueError, AttributeError) as exc:
            raise SerializationError(f"Invalid arguments to {head}") from exc
        if result.has(sp.oo, -sp.oo, sp.zoo, sp.nan):
            raise SerializationError("AST evaluates to a non-finite or undefined expression")
        return result

    return decode(payload, 0)


def _keys(node, required, optional=frozenset()):
    if not required <= node.keys() or node.keys() - required - optional:
        raise SerializationError("Unexpected or missing AST fields")


def _scalar(expr):
    if not isinstance(expr, sp.Expr) or isinstance(expr, Relational):
        raise SerializationError("Expected a scalar symbolic expression")
