"""Real Wolfram/HypExp integration tests. Never replaced by mock output.

Opt in with HYPRA_RUN_REAL=1 after configuring external dependencies.
"""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
import sympy as s

from hypra import (
    HPL,
    HPLPM,
    BackendConfig,
    EulerSum,
    HPLs6,
    NativeConfig,
    NielsenPolylog,
    backend_info,
    evaluate_hpl_boundary,
    expand_epsilon,
    hpl_at_i,
    hpl_pm_at_i,
    reduce_euler_sum,
    reduce_hpl_to_minimal_set,
    reduce_polylogs,
    transform_polylog_arguments,
)
from hypra.hpl.at_i_reduction import _PRIMITIVES
from hypra.hpl.multiple_polylog import MultiplePolylog
from hypra.hpl.special import _reduce_at_i_constants

pytestmark = pytest.mark.real_backend
E, W = s.symbols("eps w")
DOMAIN = s.And(W > 0, W < 1)


@pytest.fixture(scope="module", autouse=True)
def require_real_backend(tmp_path_factory):
    if os.environ.get("HYPRA_RUN_REAL") != "1":
        pytest.skip("Real backend is opt-in: HYPRA_RUN_REAL=1")
    with pytest.MonkeyPatch.context() as env:
        env.setenv("HYPRA_CACHE", "0")
        if not (os.environ.get("HYPRA_HYPEXP_PATH") or os.environ.get("HYPRA_HPL_PATH")):
            upstream = Path(__file__).resolve().parents[3] / "reference" / "upstream"
            hpl = upstream / "HPL-2.0" / "HPL.m"
            hypexp = upstream / "HypExp-2.0"
            if hpl.is_file() and (hypexp / "HypExp.m").is_file():
                workspace = tmp_path_factory.mktemp("upstream-hypexp") / "HypExp-2.0"
                shutil.copytree(hypexp, workspace)
                known = workspace / ("known.m" if os.name == "nt" else ".known.m")
                if not known.exists():
                    known.write_text("{KnownToOrder1}\n", encoding="utf-8")
                env.setenv("HYPRA_HPL_PATH", str(hpl))
                env.setenv("HYPRA_HYPEXP_PATH", str(workspace / "HypExp.m"))
        info = backend_info(backend="wolfram")
        assert info["available"], info["error"]
        yield


def elementary_hpl(expr):
    """Known real-interval identities for only the words used by these tests."""
    replacements = {}
    for atom in expr.atoms(HPL):
        word, z = tuple(atom.word), atom.argument
        if word == (0,):
            replacements[atom] = s.log(z)
        elif word and all(letter == 1 for letter in word):
            replacements[atom] = (-s.log(1 - z)) ** len(word) / s.factorial(len(word))
        elif word and word[-1] == 1 and all(letter == 0 for letter in word[:-1]):
            replacements[atom] = s.polylog(len(word), z)
    return expr.xreplace(replacements)


def test_p_and_q_exact_coefficients():
    p = expand_epsilon(s.hyper([1, -E], [1 - 2 * E], W), E, 2, DOMAIN, backend="wolfram")
    q = expand_epsilon(s.hyper([1, -2 * E], [1 - 2 * E], 1 - W), E, 2, DOMAIN, backend="wolfram")
    assert s.simplify(elementary_hpl(p.coefficient(1)) - s.log(1 - W)) == 0
    assert (
        s.simplify(elementary_hpl(p.coefficient(2)) + 2 * s.polylog(2, W) + s.log(1 - W) ** 2 / 2)
        == 0
    )
    assert s.simplify(elementary_hpl(q.coefficient(1)) - 2 * s.log(W)) == 0
    assert s.simplify(elementary_hpl(q.coefficient(2)) + 4 * s.polylog(2, 1 - W)) == 0
    assert p.provenance["injected_backend"] is False


def test_gamma_pole_forces_deeper_hypergeometric_order():
    q = s.hyper([1, -2 * E], [1 - 2 * E], 1 - W)
    result = expand_epsilon(s.gamma(E) * q, E, through=2, assumptions=DOMAIN, backend="wolfram")
    q_series = (
        1 + 2 * E * s.log(W) - 4 * E**2 * s.polylog(2, 1 - W) - 8 * E**3 * s.polylog(3, 1 - W)
    )
    expected = s.series(s.gamma(E) * q_series, E, 0, 3).removeO().expand()
    for power in range(-1, 3):
        assert s.simplify(elementary_hpl(result.coefficient(power)) - expected.coeff(E, power)) == 0
    assert result.provenance["plan"]["internal_order"] == 3


def test_half_integer_family_at_exact_argument():
    f = s.hyper([1, s.Rational(1, 2) + E], [s.Rational(3, 2) + E], s.Rational(1, 4))
    result = expand_epsilon(f, E, through=0, backend="wolfram")
    # 2F1(1,1/2;3/2;1/4)=2 atanh(1/2)=log(3).
    assert s.simplify(elementary_hpl(result.coefficient(0)) - s.log(3)) == 0


def test_half_integer_plus_letter_coefficient():
    half = s.S.Half
    f = s.hyper([1, half + E], [1 + half + E], W)
    upstream = expand_epsilon(f, E, through=1, assumptions=DOMAIN, backend="wolfram")
    native = expand_epsilon(f, E, through=1, assumptions=DOMAIN)
    from hypra.symbolic import reduce_known_hpl

    actual = reduce_known_hpl(upstream.coefficient(1))
    expected = native.coefficient(1)
    for point in (s.Rational(1, 4), s.Rational(4, 9)):
        difference = (actual - expected).subs(W, point)
        assert abs(difference.evalf(40)) < s.Float("1e-35")


def test_half_integer_unit_endpoint():
    f = s.hyper([s.S.Half, 1, 1 + E], [s.Rational(3, 2), 2], 1)
    upstream = expand_epsilon(f, E, through=1, backend="wolfram")
    native = expand_epsilon(f, E, through=1)
    assert s.simplify(upstream.polynomial - native.polynomial) == 0


def test_hpl_weight_four_at_i():
    expression = HPL((0, 1, 0, -1), s.I)
    upstream = expand_epsilon(expression, E, through=0, backend="wolfram")
    assert not upstream.polynomial.has(HPL)
    assert abs(upstream.polynomial.evalf(36) - expression.evalf(36)) < s.Float("1e-30")


def _hpl_reference(tmp_path, source):
    config = BackendConfig.from_env()
    info = backend_info(backend="wolfram")
    executable = config.resolved_executable()
    assert executable is not None
    script = tmp_path / "compare_hpl.wls"
    output = tmp_path / "values.json"
    script.write_text(
        r"""argv = If[$ScriptCommandLine === {}, $CommandLine, $ScriptCommandLine];
Global`$HPLPath = DirectoryName[argv[[-2]]];
Get[argv[[-2]]];
words = Flatten[Table[Tuples[{-1, 0, 1}, n], {n, 1, 6}], 1];
toPM[word_] := Replace[word, {-1 -> HPL`minus, 1 -> HPL`plus}, {1}];
fromPM[word_] := Replace[word, {HPL`minus -> -1, HPL`plus -> 1}, {1}];
"""
        + source,
        encoding="utf-8",
    )
    mode = "-file" if "wolframscript" in Path(executable).name.lower() else "-script"
    run = subprocess.run(
        [executable, mode, str(script), info["hpl_path"], str(output)],
        capture_output=True,
        text=True,
        timeout=max(config.timeout, 180),
        check=False,
    )
    assert run.returncode == 0 and output.is_file(), run.stdout + run.stderr
    return json.loads(output.read_text())


def test_all_hpl_words_through_weight_six_at_i(tmp_path):
    rows = _hpl_reference(
        tmp_path,
        r"""
pack[value_] := {ToString[Round[10^42 Re[value]], InputForm],
                 ToString[Round[10^42 Im[value]], InputForm]};
pmWords = {{-1}, {1}, {1, 1, 1, 1, 1, 1}, {-1, -1, -1, -1, -1, -1},
           {0, -1, 1, -1, 0}, {1, 0, -1, 0, 1, -1}};
words = Prepend[words, {}];
ordinary = Table[N[If[word === {}, 1, HPL`HPL[word, I]], 55], {word, words}];
pm = Table[N[HPL`HPL[toPM[word], I], 55], {word, pmWords}];
If[!AllTrue[Join[ordinary, pm], NumberQ], Exit[2]];
Export[Last[argv], {MapThread[Prepend[pack[#2], #1] &, {words, ordinary}],
                   MapThread[Prepend[pack[#2], #1] &, {pmWords, pm}]}, "RawJSON"];
Exit[0];
""",
    )
    ordinary, pm = rows
    assert len({tuple(word) for word, _, _ in ordinary}) == 1093
    assert len(pm) == 6
    constants = {}
    primitives = {HPLPM.from_indices(word, s.I) for words in _PRIMITIVES.values() for word in words}
    for items, evaluate in [
        (ordinary, lambda word: hpl_at_i(tuple(word))),
        (pm, lambda word: hpl_pm_at_i(_pm_word(word))),
    ]:
        for word, real, imag in items:
            upstream = (s.Integer(real) + s.I * s.Integer(imag)) / 10**42
            actual = evaluate(word)
            assert not actual.has(HPL, MultiplePolylog, s.Float), word
            assert actual.atoms(HPLPM) <= primitives, word
            for constant in actual.atoms(HPLPM, s.polylog):
                if constant not in constants:
                    constants[constant] = constant.evalf(55)
            numerical = actual.xreplace(constants).evalf(40)
            assert abs(numerical - upstream.evalf(40)) < s.Rational(1, 10**35), word


def _pm_word(indices):
    letters = {-1: "minus", 0: 0, 1: "plus"}
    return tuple(letters[index] for index in indices)


def test_hpl_minimal_sets_match_live_tables(tmp_path):
    rows = _hpl_reference(
        tmp_path,
        r"""
HPL`$HPLAutoConvertToKnownFunctions = False;
HPL`$HPLAutoConvertToSimplerArgument = False;
HPL`$HPLAutoLogExtract = False;
HPL`$HPLAutoProductExpand = False;
HPL`$HPLAutoReduceToMinimalSet = False;
pack[word_] := Module[{result, vars, coefficients},
  result = Expand[HPL`HPLReduceToMinimalSet[HPL`HPL[word, x]]];
  vars = DeleteDuplicates[Cases[result, h_HPL`HPL :> h, {0, Infinity}]];
  coefficients = CoefficientRules[result, vars];
  {Map[fromPM[HPL`HPLMtoA[First[#]]] &, vars],
   Map[{First[#], ToString[Last[#], InputForm]} &, coefficients]}
];
Export[Last[argv], {Table[Prepend[pack[word], word], {word, words}],
                   Table[Prepend[pack[toPM[word]], word], {word, words}]}, "RawJSON"];
Exit[0];
""",
    )
    z = s.Symbol("z")
    for items, construct in [
        (rows[0], lambda word: HPL(tuple(word), z)),
        (rows[1], lambda word: HPLPM(_pm_word(word), z)),
    ]:
        assert len({tuple(word) for word, _, _ in items}) == 1092
        for word, generators, terms in items:
            basis = [s.Dummy() for _ in generators]
            expected = s.Add(
                *(
                    s.Rational(coefficient)
                    * s.prod(atom**power for atom, power in zip(basis, powers, strict=True))
                    for powers, coefficient in terms
                )
            )
            actual = reduce_hpl_to_minimal_set(construct(word))
            assert (
                actual.xreplace(
                    {
                        construct(indices): symbol
                        for indices, symbol in zip(generators, basis, strict=True)
                    }
                ).expand()
                == expected
            ), word


@pytest.mark.parametrize(
    "upper,lower,through",
    [
        ([1 + E, 1 - 2 * E], [E], 0),
        ([2 + E, 1 - E], [3 + 2 * E], 1),
        ([E, 2 * E], [s.S.Half + 3 * E], 2),
        ([1, E, 2 * E], [s.S.Half + 3 * E, 1 + 4 * E], 2),
        ([s.S.Half + E, 2 * E], [1 + 3 * E], 1),
        ([s.S.Half + E, 2 * E, 3 * E], [1 + 4 * E, 1 + 5 * E], 2),
        ([1 + E, s.S.Half + 2 * E], [s.S.Half + 3 * E], 1),
        ([1 + E, 1 + 2 * E, s.S.Half + 3 * E], [1 + 4 * E, s.S.Half + 5 * E], 1),
        ([1 + E, s.S.Half + 2 * E, s.S.Half + 3 * E], [s.S.Half + 4 * E, s.S.Half + 5 * E], 1),
        ([s.S.Half + E, s.S.Half + 2 * E], [s.S.Half + 3 * E], 1),
        (
            [s.S.Half + E, s.S.Half + 2 * E, s.S.Half + 3 * E],
            [s.S.Half + 4 * E, s.S.Half + 5 * E],
            1,
        ),
        (
            [1 + E, 1 + 2 * E, 1 + 3 * E, s.S.Half + 4 * E],
            [1 + 5 * E, 1 + 6 * E, s.S.Half + 7 * E],
            1,
        ),
        (
            [1 + E, s.S.Half + 2 * E, s.S.Half + 3 * E, s.S.Half + 4 * E],
            [s.S.Half + 5 * E, s.S.Half + 6 * E, s.S.Half + 7 * E],
            1,
        ),
    ],
    ids=[
        "integer-pole",
        "integer-offset",
        "A-rank2",
        "A-rank3",
        "B-rank2",
        "B-rank3",
        "matched-rank2",
        "matched-one-rank3",
        "matched-rank3",
        "all-half-rank2",
        "all-half-rank3",
        "matched-one-rank4",
        "matched-three-rank4",
    ],
)
def test_hypexp_dispatch_classes_live(upper, lower, through):
    expression = s.hyper(upper, lower, W, evaluate=False)
    native = expand_epsilon(expression, E, through, DOMAIN, config=NativeConfig(use_cache=False))
    upstream = expand_epsilon(expression, E, through, DOMAIN, backend="wolfram")
    first = -native.provenance["plan"]["pole_bound"]
    for power in range(first, through + 1):
        difference = native.coefficient(power) - upstream.coefficient(power)
        for point in (s.Rational(1, 4), s.Rational(4, 9)):
            assert abs(difference.subs(W, point).evalf(40)) < s.Rational(1, 10**32)


def test_polylog_argument_relations_live(tmp_path):
    rows = _hpl_reference(
        tmp_path,
        r"""
indices = {{1,1},{2,1},{3,1},{4,1},{1,2},{2,2},{1,3},{2,3},{3,2}};
z = 2/5;
arguments = {z,1-z,z/(z-1),1-1/z,1/z,1/(1-z)};
rows = Flatten[Table[
  value = N[If[np[[2]] == 1, PolyLog[np[[1]]+1, arguments[[j]]],
                              PolyLog[np[[1]],np[[2]],arguments[[j]]]], 48];
  If[!NumberQ[value], Exit[2]];
  {np[[1]],np[[2]],j,ToString[Round[10^40 Re[value]], InputForm],
                       ToString[Round[10^40 Im[value]], InputForm]},
  {np,indices},{j,Length[arguments]}], 1];
Export[Last[argv], rows, "RawJSON"];
Exit[0];
""",
    )
    assert len(rows) == 54
    arguments = [W, 1 - W, W / (W - 1), 1 - 1 / W, 1 / W, 1 / (1 - W)]
    for n, p, index, real, imag in rows:
        source = (
            s.polylog(n + 1, arguments[index - 1])
            if p == 1
            else NielsenPolylog(n, p, arguments[index - 1])
        )
        actual = transform_polylog_arguments(source, W, domain="0<z<1", boundary="lower")
        expected = (s.Integer(real) + s.I * s.Integer(imag)) / 10**40
        assert abs((actual.subs(W, s.Rational(2, 5)) - expected).evalf(36)) < s.Rational(1, 10**30)
        named = reduce_polylogs(source, W, domain="0<z<1", boundary="lower")
        assert abs((named.subs(W, s.Rational(2, 5)) - expected).evalf(36)) < s.Rational(1, 10**30)


def test_mixed_weight_seven_real_boundaries_live(tmp_path):
    rows = _hpl_reference(
        tmp_path,
        r"""
word = {1,0,-1,0,1,0,-1};
rows = Flatten[Table[
  value = N[HPL`HPL[word,x,HPL`AnalyticContinuationSign->side],48];
  If[!NumberQ[value], Exit[2]];
  {x,side,ToString[Round[10^40 Re[value]],InputForm],
          ToString[Round[10^40 Im[value]],InputForm]},
  {x,{-2,2}},{side,{-1,1}}],1];
Export[Last[argv], rows, "RawJSON"];
Exit[0];
""",
    )
    assert len(rows) == 4
    for point, side, real, imag in rows:
        result = evaluate_hpl_boundary(
            (1, 0, -1, 0, 1, 0, -1),
            point,
            region="x>1" if point > 0 else "x<-1",
            boundary="upper" if side > 0 else "lower",
            dps=35,
        )
        expected = (s.Integer(real) + s.I * s.Integer(imag)) / 10**40
        actual = s.Float(str(result.value.real), 40) + s.I * s.Float(str(result.value.imag), 40)
        assert abs((actual - expected).evalf(35)) < s.Rational(1, 10**30)


def test_all_colored_endpoint_constants_through_weight_six_live(tmp_path):
    rows = _hpl_reference(
        tmp_path,
        r"""
compositions[0] = {{}};
compositions[n_] := compositions[n] = Join @@ Table[
  Prepend[#, k]& /@ compositions[n-k], {k,1,n}];
rows = {};
Do[
  Do[
    Do[
      If[First[indices] > 1 || First[colors] == -1,
        (* HPL 2's generic one-index MZV rule is indeterminate at {-1};
           use the equivalent convergent HPL endpoint for that one value. *)
        value = N[If[indices == {1}, -HPL`HPL[{-1},1],
                     HPL`MZV[indices*colors]],48];
        If[!NumberQ[value], Print[{indices,colors,value}]; Exit[2]];
        AppendTo[rows,{indices,colors,ToString[Round[10^40 value],InputForm]}]],
      {colors,Tuples[{-1,1},Length[indices]]}],
    {indices,compositions[weight]}],
  {weight,1,6}];
Export[Last[argv],rows,"RawJSON"];
Exit[0];
""",
    )
    assert len(rows) == 485
    evaluated = {}
    for indices, colors, value in rows:
        expression = reduce_euler_sum(EulerSum(tuple(indices), tuple(colors)))
        for constant in expression.atoms(EulerSum, s.polylog):
            if constant not in evaluated:
                evaluated[constant] = constant.evalf(55)
        actual = expression.xreplace(evaluated).evalf(40)
        expected = s.Integer(value) / 10**40
        assert abs(actual - expected) < s.Rational(1, 10**32), (indices, colors)


def _decode_at_i_constant(tree):
    """Decode a restricted exact constant AST, never source-language code."""
    kind, *arguments = tree
    if kind == "integer":
        return s.Integer(arguments[0])
    if kind == "rational":
        return s.Rational(*arguments)
    if kind == "symbol":
        return {"Pi": s.pi, "Catalan": s.Catalan, "HPLs6": HPLs6}[arguments[0]]
    if kind == "hpli":
        word = tuple(arguments[0])
        return (
            hpl_pm_at_i(word) if any(isinstance(letter, str) for letter in word) else hpl_at_i(word)
        )
    decoded = [_decode_at_i_constant(argument) for argument in arguments]
    if kind == "complex":
        return decoded[0] + s.I * decoded[1]
    if kind == "Zeta" and len(decoded) == 2:
        order, offset = decoded
        if (
            isinstance(order, s.Integer)
            and order >= 2
            and offset in (s.Rational(1, 4), s.Rational(3, 4))
        ):
            # Separate odd residue classes: A+B=(4^n-2^n)*zeta(n),
            # A-B=4^n*beta(n). This changes only the constant coordinates.
            total = (4**order - 2**order) * s.zeta(order)
            if order == 2:
                quarter = s.pi**2 + 8 * s.Catalan
            elif order % 2:
                half = (order - 1) // 2
                beta = (
                    (-1) ** half
                    * s.euler(2 * half)
                    * s.pi**order
                    / (4 ** (half + 1) * s.factorial(2 * half))
                )
                quarter = (total + 4**order * beta) / 2
            else:
                quarter = s.polygamma(order - 1, s.Rational(1, 4)) / s.factorial(order - 1)
            return quarter if offset == s.Rational(1, 4) else total - quarter
    if kind == "PolyGamma" and len(decoded) == 2:
        order, offset = decoded
        if isinstance(order, s.Integer) and order >= 1 and offset == s.Rational(3, 4):
            # The same residue-class sum in polygamma coordinates.
            total = (
                (-1) ** (order + 1)
                * s.factorial(order)
                * (4 ** (order + 1) - 2 ** (order + 1))
                * s.zeta(order + 1)
            )
            return total - s.polygamma(order, s.Rational(1, 4))
    functions = {
        "Plus": s.Add,
        "Times": s.Mul,
        "Power": s.Pow,
        "Log": s.log,
        "PolyLog": s.polylog,
        "PolyGamma": s.polygamma,
        "Zeta": s.zeta,
    }
    return functions[kind](*decoded)


def test_at_i_named_constant_reductions_match_live_tables_exactly(tmp_path):
    rows = _hpl_reference(
        tmp_path,
        r"""
encode[x_Integer] := {"integer",ToString[x,InputForm]};
encode[x_Rational] := {"rational",ToString[Numerator[x],InputForm],
                                  ToString[Denominator[x],InputForm]};
encode[x_Complex] := {"complex",encode[Re[x]],encode[Im[x]]};
encode[x_Symbol] := If[MemberQ[{"Pi","Catalan","HPLs6"},SymbolName[x]],
                      {"symbol",SymbolName[x]},Print[x];Exit[2]];
encode[HPL`HPLI[word_]] := {"hpli",word /. {HPL`plus->"plus",HPL`minus->"minus"}};
encode[x_] := Module[{head=SymbolName[Head[x]]},
  If[!MemberQ[{"Plus","Times","Power","Log","PolyLog","PolyGamma","Zeta"},head],
     Print[InputForm[x]];Exit[2]];
  Prepend[encode /@ (List @@ x),head]];
selected = Select[words,Length[#]>=4&];
rows = Table[{word,encode[HPL`HPL[toPM[word],I]]},{word,selected}];
Export[Last[argv],rows,"RawJSON"];
Exit[0];
""",
    )
    assert len(rows) == 1053
    for word, tree in rows:
        expected = _reduce_at_i_constants(_decode_at_i_constant(tree))
        actual = hpl_pm_at_i(_pm_word(word))
        assert s.expand(actual - expected) == 0, word


@pytest.mark.parametrize("weight", [4, 5, 6])
def test_at_i_primitive_round_trip_through_real_bridge(weight):
    primitive = HPLPM.from_indices(_PRIMITIVES[weight][0], s.I)
    expression = (1 + E) * primitive
    native = expand_epsilon(expression, E, 1, config=NativeConfig(use_cache=False))
    upstream = expand_epsilon(expression, E, 1, backend="wolfram")
    assert upstream.polynomial == native.polynomial
    assert upstream.coefficient(0) == upstream.coefficient(1) == primitive
