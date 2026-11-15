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

from hypra import HPL, BackendConfig, backend_info, expand_epsilon, hpl_at_i

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


def test_all_hpl_words_through_weight_four_at_i(tmp_path):
    config = BackendConfig.from_env()
    info = backend_info(backend="wolfram")
    executable = config.resolved_executable()
    assert executable is not None
    script = tmp_path / "compare_i.wls"
    output = tmp_path / "values.json"
    script.write_text(
        r"""argv = If[$ScriptCommandLine === {}, $CommandLine, $ScriptCommandLine];
Global`$HPLPath = DirectoryName[argv[[-2]]];
Get[argv[[-2]]];
words = Join[{{}}, Flatten[Table[Tuples[{-1, 0, 1}, n], {n, 1, 4}], 1]];
values = Table[N[If[word === {}, 1, HPL`HPL[word, I]], 55], {word, words}];
If[!AllTrue[values, NumberQ], Exit[2]];
rows = MapThread[{#1, ToString[Round[10^42 Re[#2]], InputForm],
                      ToString[Round[10^42 Im[#2]], InputForm]} &, {words, values}];
Export[Last[argv], rows, "RawJSON"];
Exit[0];
""",
        encoding="utf-8",
    )
    mode = "-file" if "wolframscript" in Path(executable).name.lower() else "-script"
    run = subprocess.run(
        [executable, mode, str(script), info["hpl_path"], str(output)],
        capture_output=True,
        text=True,
        timeout=config.timeout,
        check=False,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    rows = json.loads(output.read_text())
    assert len({tuple(word) for word, _, _ in rows}) == 121
    for word, real, imag in rows:
        upstream = (s.Integer(real) + s.I * s.Integer(imag)) / 10**42
        actual = hpl_at_i(tuple(word))
        assert not actual.has(HPL, s.Float), word
        assert abs(actual.evalf(40) - upstream.evalf(40)) < s.Float("1e-35"), word
