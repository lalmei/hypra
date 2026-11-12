"""Convergent endpoint identities, exact MZV transport and divergence guards."""

import json

import mpmath as mp
import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl_algebra import shuffle_words
from hypra.hpl_endpoints import dual_indices, hpl_at_one, reduce_mzv
from hypra.serialization import SerializationError, decode_expr, encode_expr
from hypra.symbolic import MZV


@pytest.mark.parametrize("indices", [(2,), (2, 1), (3, 2), (2, 1, 1, 1)])
def test_mzv_lossless_json_and_latex(indices):
    value = MZV(indices)
    assert decode_expr(json.loads(json.dumps(encode_expr(value)))) == value
    assert value.indices == s.Tuple(*indices)
    assert s.diff(value, s.Symbol("z")) == 0
    assert s.latex(value) == r"\zeta\!\left(" + ",".join(map(str, indices)) + r"\right)"
    assert value.is_real is True
    assert value.is_positive is True


def test_mzv_empty_and_symbolic_combination():
    assert MZV(()) == 1
    value = s.Rational(2, 3) * MZV((3, 2)) + s.zeta(3) * MZV((2, 1))
    assert decode_expr(encode_expr(value)) == value


@pytest.mark.parametrize(
    "indices", [(1,), (1, 2), (0,), (-2,), (2, 0), (s.Rational(3, 2),), (2.0,)]
)
def test_nonconvergent_or_noninteger_mzv_is_rejected(indices):
    with pytest.raises((TypeError, ValueError)):
        MZV(indices)


@pytest.mark.parametrize(
    "ast",
    [
        {"head": "MZV", "args": []},
        {"head": "MZV", "args": [{"head": "Integer", "value": "2"}]},
        {"head": "MZV", "args": [{"head": "List", "args": [{"head": "Integer", "value": "1"}]}]},
    ],
)
def test_invalid_mzv_wire_data_is_rejected(ast):
    with pytest.raises(SerializationError):
        decode_expr(ast)


def test_admissible_hpl_to_mzv_is_lossless_before_reduction():
    assert hpl_at_one((0, 1), reduced=False) == MZV((2,))
    assert hpl_at_one((0, 0, 1, 0, 1), reduced=False) == MZV((3, 2))
    assert hpl_at_one((0, 1, 1), reduced=False) == MZV((2, 1))
    assert hpl_at_one(()) == 1
    assert hpl_at_one((0, 0, 0)) == 0


def test_duality_and_derived_low_weight_reductions():
    assert dual_indices((2, 1, 1)) == (4,)
    assert dual_indices((3, 2)) == (2, 2, 1)
    assert dual_indices(dual_indices((3, 2))) == (3, 2)
    assert reduce_mzv(MZV((2, 1))) == s.zeta(3)
    assert reduce_mzv(MZV((2, 1, 1))) == s.zeta(4)
    assert reduce_mzv(MZV((3, 1))) == s.zeta(4) / 4
    assert reduce_mzv(MZV((2, 2))) == 3 * s.zeta(4) / 4
    assert reduce_mzv(MZV((3, 2))) == MZV((3, 2))
    assert reduce_mzv(MZV((25, 2))) == MZV((25, 2))
    assert reduce_mzv(MZV((25,))) == s.zeta(25)


def test_weight_four_shuffle_stuffle_consistency():
    shuffled = shuffle_words((0, 1), (0, 1))
    endpoint_product = sum(count * hpl_at_one(word) for word, count in shuffled.items())
    assert s.simplify(endpoint_product - s.zeta(2) ** 2) == 0
    assert s.simplify(2 * reduce_mzv(MZV((2, 2))) + s.zeta(4) - s.zeta(2) ** 2) == 0


@pytest.mark.parametrize("zeros", range(1, 7))
def test_convergent_leading_one_with_zero_tail(zeros):
    # Direct integral log(t)^m/(m!*(1-t)) = (-1)^m zeta(m+1).
    assert s.simplify(hpl_at_one((1,) + (0,) * zeros) - (-1) ** zeros * s.zeta(zeros + 1)) == 0


def test_trailing_zero_shuffle_preserves_complete_limit():
    assert hpl_at_one((0, 1, 0), reduced=False) == -2 * MZV((3,))
    assert hpl_at_one((0, 1, 0, 0)) == 3 * s.zeta(4)


@pytest.mark.parametrize("word", [(1,), (1, 1), (1, 0, 1), (1, 1, 0), (1, 0, 1, 0)])
def test_divergent_endpoint_not_regularized(word):
    with pytest.raises(UnsupportedExpressionError, match="Divergent"):
        hpl_at_one(word)


def test_other_alphabets_and_bad_input_fail_explicitly():
    with pytest.raises(UnsupportedExpressionError, match="alphabet"):
        hpl_at_one((0, -1))
    with pytest.raises(TypeError):
        hpl_at_one("0,1")


def test_secondary_endpoint_quadrature_agrees_with_exact_identities():
    # Independent scalar integrals, not numerical fitting or implementation data.
    with mp.workdps(35):
        weight3 = mp.quad(lambda t: mp.log1p(-t) ** 2 / (2 * t), [0, mp.mpf(".5"), 1])
        weight4 = mp.quad(lambda t: -mp.log(t) * mp.log1p(-t) ** 2 / (2 * t), [0, mp.mpf(".5"), 1])
        assert abs(weight3 - mp.zeta(3)) < mp.mpf("1e-30")
        assert abs(weight4 - mp.zeta(4) / 4) < mp.mpf("1e-30")


def test_exact_index_validation_cannot_be_bypassed_by_function_cache():
    from hypra.symbolic import HPL

    z = s.Symbol("z")
    MZV((2,))
    HPL((1,), z)
    with pytest.raises(ValueError):
        MZV((2.0,))
    with pytest.raises(ValueError):
        HPL((1.0,), z)
    with pytest.raises(ValueError):
        MZV((1,), evaluate=False)


def test_endpoint_public_exports():
    from hypra import MZV as PublicMZV
    from hypra import hpl_at_one as public_endpoint

    assert PublicMZV((3, 2)) == MZV((3, 2))
    assert public_endpoint((0, 1)) == s.zeta(2)
