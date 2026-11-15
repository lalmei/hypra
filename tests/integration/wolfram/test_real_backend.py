"""Real Wolfram/HypExp integration tests. Never replaced by mock output.

Opt in with HYPRA_RUN_REAL=1 after configuring external dependencies.
"""

import os

import pytest
import sympy as s

from hypra import HPL, backend_info, expand_epsilon

pytestmark = pytest.mark.real_backend
E, W = s.symbols("eps w")
DOMAIN = s.And(W > 0, W < 1)


@pytest.fixture(scope="module", autouse=True)
def require_real_backend():
    if os.environ.get("HYPRA_RUN_REAL") != "1":
        pytest.skip("Real backend is opt-in: HYPRA_RUN_REAL=1")
    info = backend_info(backend="wolfram")
    assert info["available"], info["error"]


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
