"""Persistent templates must preserve the native expansion contract."""

import sqlite3

import pytest
import sympy as s

from hypra import ExpansionLibrary, NativeConfig, UnsupportedExpressionError, expand_epsilon

E, Z = s.symbols("eps z")
DOMAIN = s.And(Z > 0, Z < 1)


@pytest.fixture(scope="module")
def library(tmp_path_factory):
    value = ExpansionLibrary(tmp_path_factory.mktemp("expansion-library") / "families.sqlite")
    value.add_to_lib([1, 0], [1], through=3)
    return value


def test_missing_library_query_does_not_create_file(tmp_path):
    path = tmp_path / "absent" / "families.sqlite"
    assert not ExpansionLibrary(path).is_known_to_order([0, 1], [1], through=2)
    assert not path.parent.exists()


def test_base_order_and_registered_depth(library):
    reopened = ExpansionLibrary(library.path)
    assert reopened.is_known_to_order([0, 1], [1], through=3)
    assert reopened.is_known_to_order([1, 0], [1], through=1)
    assert not reopened.is_known_to_order([0, 1], [1], through=4)
    assert not reopened.is_known_to_order([0, 2], [1], through=1)


@pytest.mark.parametrize("slopes", [(1, 2, 3), (0, -2, 0), (1, 1, 1)])
def test_specialized_template_matches_native_kernel(library, slopes):
    a, b, c = slopes
    atom = s.hyper([a * E, 1 + b * E], [1 + c * E], Z, evaluate=False)
    expected = expand_epsilon(atom, E, 2, DOMAIN, config=NativeConfig(use_cache=False))
    actual = expand_epsilon(
        atom,
        E,
        2,
        DOMAIN,
        config=NativeConfig(use_cache=False, expansion_library=str(library.path)),
    )
    assert s.simplify(actual.polynomial - expected.polynomial) == 0


def test_template_reuse_respects_whole_expression_pole_depth(library, monkeypatch):
    from hypra.native import backend

    atom = s.hyper([1, -E], [1 - 2 * E], Z)
    expression = s.gamma(E) * atom
    expected = expand_epsilon(expression, E, 2, DOMAIN, config=NativeConfig(use_cache=False))

    def no_kernel(*args, **kwargs):
        raise AssertionError("Registered template should be used through internal order three")

    monkeypatch.setattr(backend, "expand_atom", no_kernel)
    actual = expand_epsilon(
        expression,
        E,
        2,
        DOMAIN,
        config=NativeConfig(use_cache=False, expansion_library=str(library.path)),
    )
    assert actual.provenance["plan"]["internal_order"] == 3
    assert s.simplify(actual.polynomial - expected.polynomial) == 0


@pytest.mark.parametrize("parameter", [s.Rational(1, 3), s.sqrt(2)])
def test_out_of_template_scope_keeps_native_elementary_expansion(tmp_path, parameter):
    atom = s.hyper([parameter], [], s.S.Half)
    expected = expand_epsilon(atom, E, 0, config=NativeConfig(use_cache=False))
    path = tmp_path / "missing.sqlite"
    actual = expand_epsilon(
        atom, E, 0, config=NativeConfig(use_cache=False, expansion_library=str(path))
    )
    assert actual.polynomial == expected.polynomial
    assert not path.exists()


def test_half_integer_library_miss_uses_native_kernel(library):
    atom = s.hyper([1, s.S.Half + E], [s.Rational(3, 2) + E], s.Rational(1, 4))
    expected = expand_epsilon(atom, E, 1, config=NativeConfig(use_cache=False))
    actual = expand_epsilon(
        atom,
        E,
        1,
        config=NativeConfig(use_cache=False, expansion_library=str(library.path)),
    )
    assert s.simplify(actual.polynomial - expected.polynomial) == 0


def test_failed_upgrade_preserves_existing_template(library, monkeypatch):
    from hypra import api

    with sqlite3.connect(library.path) as database:
        before = database.execute("SELECT * FROM expansions").fetchall()

    def fail(*args, **kwargs):
        raise UnsupportedExpressionError("Deliberate algebra budget failure")

    monkeypatch.setattr(api, "expand_epsilon", fail)
    with pytest.raises(UnsupportedExpressionError, match="budget"):
        library.add_to_lib([0, 1], [1], through=4)
    with sqlite3.connect(library.path) as database:
        assert database.execute("SELECT * FROM expansions").fetchall() == before
