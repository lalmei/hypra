import json

import pytest
import sympy as sp

from hypra.serialization import SerializationError, decode_expr, encode_expr
from hypra.symbolic import HPL

EPS, Z = sp.symbols("eps z")


@pytest.mark.parametrize("expr", [
    sp.Integer(-12345678901234567890), sp.Rational(-37, 19),
    sp.pi, sp.E, sp.I, sp.EulerGamma, sp.Catalan, sp.true, sp.false,
    EPS, sp.Symbol("positive_z", positive=True),
    2 * EPS**-2 + sp.gamma(1 - EPS) * sp.hyper((1, -EPS), (1 - 2*EPS,), Z),
    sp.log(Z) + sp.polylog(2, Z) + sp.zeta(3) + sp.exp(Z),
    sp.polygamma(2, Z), sp.zeta(2, Z),
    HPL((0, 1, -1), Z), HPL((), Z), HPL((0, 0), 1 - Z),
    sp.And(Z > 0, Z < 1), sp.Or(Z < 0, Z > 1), sp.Not(sp.Eq(Z, 0)),
    sp.Eq(Z, 1, evaluate=False), sp.Ne(Z, 0, evaluate=False),
    sp.Le(Z, 1, evaluate=False), sp.Ge(Z, 0, evaluate=False),
])
def test_exact_json_roundtrip(expr):
    ast = json.loads(json.dumps(encode_expr(expr)))
    assert decode_expr(ast) == expr


def test_preserves_uncancelled_hypergeometric_parameters():
    expr = sp.hyper((1, EPS), (EPS,), Z, evaluate=False)
    result = decode_expr(encode_expr(expr))
    assert result.ap == (1, EPS)
    assert result.bq == (EPS,)


@pytest.mark.parametrize("expr", [
    0.5, "1+z", sp.Float("0.1"), sp.oo, -sp.oo, sp.zoo, sp.nan,
    sp.Integral(Z, Z), sp.Derivative(sp.Function("f")(Z), Z),
    sp.Function("unknown")(Z), sp.Order(EPS**3),
    sp.Symbol("System`Run"), sp.Symbol("x;Run[command]"),
    sp.Symbol("epsilon-with-hyphens"), sp.Symbol("nc", commutative=False),
    sp.Dummy("x"), sp.Wild("x"),
])
def test_unsupported_input_fails_closed(expr):
    with pytest.raises(SerializationError):
        encode_expr(expr)


@pytest.mark.parametrize("ast", [
    None, [], {}, {"head": 1}, {"head": "Run", "args": []},
    {"head": "Integer", "value": 1},
    {"head": "Integer", "value": "01"},
    {"head": "Integer", "value": "1;Quit[]"},
    {"head": "Integer", "value": "1", "code": "Run[]"},
    {"head": "Rational", "numerator": "1", "denominator": "0"},
    {"head": "Rational", "numerator": "1", "denominator": "-2"},
    {"head": "Symbol", "name": "x", "assumptions": {"positive": "true"}},
    {"head": "Symbol", "name": "x", "assumptions": {"malicious": True}},
    {"head": "Symbol", "name": "x", "assumptions": {"positive": True, "negative": True}},
    {"head": "Symbol", "name": "Global`x"},
    {"head": "Symbol", "name": "nc", "assumptions": {"commutative": False}},
    {"head": "Pi", "args": []},
    {"head": "Pow", "args": []},
    {"head": "Gamma", "args": "x"},
    {"head": "Gamma", "args": [{"head": "List", "args": []}]},
    {"head": "Gamma", "args": [{"head": "Integer", "value": "0"}]},
    {"head": "And", "args": [{"head": "Integer", "value": "2"}]},
    {"head": "HypergeometricPFQ", "args": [
        {"head": "Integer", "value": "1"},
        {"head": "List", "args": []}, {"head": "Symbol", "name": "z"},
    ]},
    {"head": "HPL", "args": [
        {"head": "List", "args": [{"head": "Integer", "value": "2"}]},
        {"head": "Symbol", "name": "z"},
    ]},
])
def test_malformed_or_unsafe_ast_fails_closed(ast):
    with pytest.raises(SerializationError):
        decode_expr(ast)


def test_deep_ast_rejected():
    ast = {"head": "Symbol", "name": "z"}
    for _ in range(130):
        ast = {"head": "Exp", "args": [ast]}
    with pytest.raises(SerializationError, match="depth"):
        decode_expr(ast)


def test_hpl_compressed_convention_and_latex():
    expr = HPL.from_compressed((2, -2, 0), Z)
    assert expr.word == sp.Tuple(0, 1, 0, -1, 0)
    assert expr.argument == Z
    assert sp.latex(expr) == r"H_{0,1,0,-1,0}\!\left(z\right)"
    assert HPL.from_compressed((1,), Z) != sp.log(1-Z)


@pytest.mark.parametrize("word", [(2,), (sp.Rational(1, 2),), (sp.Symbol("a"),)])
def test_hpl_rejects_ambiguous_expanded_word(word):
    with pytest.raises(ValueError):
        HPL(word, Z)


def test_compressed_hpl_rejects_noninteger_and_excessive_word():
    with pytest.raises(TypeError):
        HPL.from_compressed((1.0,), Z)
    with pytest.raises(ValueError):
        HPL.from_compressed((10001,), Z)
