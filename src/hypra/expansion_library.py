"""Persistent native expansions indexed by integer hypergeometric bases."""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterable
from contextlib import closing
from pathlib import Path
from typing import cast

import sympy as s

from .errors import IncompleteExpansionError, UnsupportedExpressionError
from .serialization import decode_expr, encode_expr


def _partition(bases: Iterable[int | s.Integer]) -> tuple[int, ...]:
    values = tuple(bases)
    if any(type(value) is not int and not isinstance(value, s.Integer) for value in values):
        raise TypeError("Expansion library bases must be exact integers")
    return tuple(sorted(map(int, values)))


def _family(upper, lower):
    a, b = _partition(upper), _partition(lower)
    if len(a) != len(b) + 1:
        raise ValueError("Expansion library requires p upper and p-1 lower bases")
    return a, b


def _check_order(through: int) -> None:
    if type(through) is not int or through < 0:
        raise ValueError("Expansion library order must be a nonnegative integer")


def _template_symbols(a, b):
    return (
        s.Symbol("libraryArgument", positive=True),
        s.Symbol("libraryEpsilon"),
        tuple(s.Symbol(f"libraryUpper{i}") for i in range(len(a))),
        tuple(
            s.Symbol(f"libraryLower{i}", nonzero=True)
            if base <= 0
            else s.Symbol(f"libraryLower{i}")
            for i, base in enumerate(b)
        ),
    )


class ExpansionLibrary:
    """An on-disk library of exact symbolic regulator expansions.

    Upper and lower bases are sorted separately. Each entry records the largest
    inclusive order registered for the family, with independent symbolic slopes.
    Pass the same SQLite file to NativeConfig(expansion_library=...) to reuse it.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path).expanduser()

    def _entry(self, key: str, through: int, *, polynomial: bool = False):
        if not self.path.exists():
            return None
        fields = "polynomial" if polynomial else "1"
        with closing(sqlite3.connect(self.path.resolve().as_uri() + "?mode=ro", uri=True)) as db:
            if (
                db.execute(
                    "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'expansions'"
                ).fetchone()
                is None
            ):
                return None
            return db.execute(
                f"SELECT {fields} FROM expansions WHERE family = ? AND known_order >= ?",
                (key, through),
            ).fetchone()

    def is_known_to_order(
        self,
        upper: Iterable[int | s.Integer],
        lower: Iterable[int | s.Integer],
        through: int,
    ) -> bool:
        """Return whether this family is registered through the requested order."""
        _check_order(through)
        return self._entry(json.dumps(_family(upper, lower)), through) is not None

    def add_to_lib(
        self,
        upper: Iterable[int | s.Integer],
        lower: Iterable[int | s.Integer],
        through: int,
    ) -> None:
        """Compute and persist a family's expansion if the requested order is new."""
        from .api import expand_epsilon
        from .native import NativeConfig

        _check_order(through)
        a, b = _family(upper, lower)
        key = json.dumps((a, b))
        if self._entry(key, through) is not None:
            return
        if any(abs(base) >= 10 for base in (*a, *b)):
            raise ValueError("Expansion library requires absolute base values below ten")
        z, epsilon, upper_slopes, lower_slopes = _template_symbols(a, b)
        h = s.hyper(
            [base + slope * epsilon for base, slope in zip(a, upper_slopes, strict=True)],
            [base + slope * epsilon for base, slope in zip(b, lower_slopes, strict=True)],
            z,
            evaluate=False,
        )
        expansion = expand_epsilon(h, epsilon, through, z < 1, config=NativeConfig(use_cache=False))
        payload = json.dumps(encode_expr(expansion.polynomial))
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with closing(sqlite3.connect(self.path)) as db, db:
            db.execute(
                """CREATE TABLE IF NOT EXISTS expansions (
                    family TEXT PRIMARY KEY, known_order INTEGER NOT NULL, polynomial TEXT NOT NULL
                )"""
            )
            db.execute(
                """INSERT INTO expansions VALUES (?, ?, ?)
                   ON CONFLICT(family) DO UPDATE SET
                       known_order = excluded.known_order, polynomial = excluded.polynomial
                   WHERE excluded.known_order > expansions.known_order""",
                (key, through, payload),
            )

    def _lookup_atom(self, h: s.hyper, epsilon: s.Symbol, through: int) -> s.Expr | None:
        from .native.planning import _affine_rational, laurent_coefficient

        # Templates are an optional optimization. Requests outside their domain
        # must still reach the native dispatcher, which also handles elementary
        # regulator-independent functions with other rational bases.
        if len(h.ap) != len(h.bq) + 1:
            return None
        try:
            upper = sorted((_affine_rational(p, epsilon) for p in h.ap), key=lambda pair: pair[0])
            lower = sorted((_affine_rational(p, epsilon) for p in h.bq), key=lambda pair: pair[0])
        except UnsupportedExpressionError:
            return None
        if any(not base.is_Integer for base, _ in (*upper, *lower)):
            return None
        a, b = _family((base for base, _ in upper), (base for base, _ in lower))
        row = self._entry(json.dumps((a, b)), through, polynomial=True)
        if row is None:
            return None
        stored = decode_expr(json.loads(row[0]))
        if not isinstance(stored, s.Expr) or stored.has(s.hyper):
            raise IncompleteExpansionError("Invalid stored hypergeometric expansion")
        z, regulator, upper_slopes, lower_slopes = _template_symbols(a, b)
        poles = sum(base <= 0 for base in b)
        polynomial = s.Add(
            *(
                laurent_coefficient(stored, regulator, power) * regulator**power
                for power in range(-poles, through + 1)
            )
        )
        replacements = dict(
            zip(upper_slopes + lower_slopes, (slope for _, slope in (*upper, *lower)), strict=True)
        )
        return cast(s.Expr, polynomial.xreplace(replacements | {z: h.argument, regulator: epsilon}))
