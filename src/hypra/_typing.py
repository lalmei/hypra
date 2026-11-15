"""Input types shared by the public symbolic and numerical APIs."""

from collections.abc import Iterable
from typing import Literal, cast

import sympy as s

type ExactInteger = int | s.Integer
type IntegerWord = (
    tuple[ExactInteger, ...] | list[int] | list[s.Integer] | list[ExactInteger] | s.Tuple
)
type ExactScalar = s.Expr | int
type Interval = Literal["0<z<1"]
type Boundary = Literal["upper", "lower"]
type RealRegion = Literal["0<x<1", "-1<x<0", "x>1", "x<-1"]


def as_integer_tuple(word: IntegerWord) -> tuple[int, ...]:
    """Convert a word whose entries have been validated as exact integers."""
    # SymPy Tuple stores its entries as Basic, regardless of their actual type.
    return tuple(map(int, cast(Iterable[ExactInteger], word)))
