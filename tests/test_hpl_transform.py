import itertools

import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import hpl_series, shuffle_expand
from hypra.hpl_transform import negate_hpl, reflect_hpl
from hypra.symbolic import HPL, reduce_known_hpl

Z = s.Symbol("z", positive=True)
DOMAIN = {"domain": "0<z<1"}


@pytest.mark.parametrize(
    "word", [w for n in range(1, 4) for w in itertools.product((0, 1), repeat=n)]
)
def test_reflection_differential_and_involution(word):
    result = reflect_hpl(word, Z, **DOMAIN)
    tail = reflect_hpl(word[1:], Z, **DOMAIN)
    kernel = -1 / (1 - Z) if word[0] == 0 else -1 / Z
    residual = reduce_known_hpl(shuffle_expand(s.diff(result, Z) - kernel * tail, Z))
    assert s.simplify(s.expand_log(residual, force=True)) == 0
    # Reflect the returned expression again, including elementary logarithms.
    # Only original log(z)
    # occurrences need substitution, not any logs introduced by replacements.
    dummy = s.Dummy("original_log")
    original = result.xreplace({s.log(Z): dummy})
    double = original.xreplace({h: reflect_hpl(h.word, Z, **DOMAIN) for h in original.atoms(HPL)})
    double = double.subs(dummy, s.log(1 - Z))
    assert (
        s.simplify(
            s.expand_log(reduce_known_hpl(shuffle_expand(double - HPL(word, Z), Z)), force=True)
        )
        == 0
    )


def test_dilog_reflection_and_weight_three_endpoint():
    residual = reflect_hpl((0, 1), Z, **DOMAIN) - (
        s.zeta(2) - s.polylog(2, Z) - s.log(Z) * s.log(1 - Z)
    )
    assert s.simplify(shuffle_expand(residual, Z)) == 0
    assert s.limit(hpl_series(reflect_hpl((0, 1, 1), Z, **DOMAIN), Z, 2), Z, 0) == s.zeta(3)


@pytest.mark.parametrize("word", [(1,), (-1,), (0, 1), (1, -1), (0, -1, 1)])
def test_negation_matches_exact_origin_series(word):
    actual = hpl_series(negate_hpl(word, Z, **DOMAIN), Z, 6)
    reference = hpl_series(HPL(word, Z), Z, 6).subs(Z, -Z)
    assert s.expand(actual - reference) == 0


def test_explicit_branch_and_fail_closed():
    assert negate_hpl((0,), Z, **DOMAIN) == s.log(Z) + s.I * s.pi
    for function, word in [(reflect_hpl, (-1,)), (reflect_hpl, (0,) * 8 + (1,))]:
        with pytest.raises(UnsupportedExpressionError):
            function(word, Z, **DOMAIN)
    with pytest.raises(UnsupportedExpressionError):
        reflect_hpl((1,), Z)


def test_mixed_zero_negation_has_correct_log_monodromy():
    result = negate_hpl((1, 0), Z, **DOMAIN)
    expected = -(s.log(Z) + s.I * s.pi) * HPL((-1,), Z) + HPL((0, -1), Z)
    assert s.expand(result - expected) == 0


@pytest.mark.parametrize("word", [(0, 0, 1, 1), (0, 1, 0, 1, 1), (1, 0, 0, 1)])
def test_higher_weight_reflection_derivative(word):
    result = reflect_hpl(word, Z, **DOMAIN)
    kernel = -1 / (1 - Z) if word[0] == 0 else -1 / Z
    residual = shuffle_expand(s.diff(result, Z) - kernel * reflect_hpl(word[1:], Z, **DOMAIN), Z)
    assert s.simplify(residual) == 0
