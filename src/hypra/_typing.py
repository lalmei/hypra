"""Input types shared by the public symbolic and numerical APIs."""

from typing import Literal

import sympy as s

type ExactInteger = int | s.Integer
type IntegerWord = (
    tuple[ExactInteger, ...] | list[int] | list[s.Integer] | list[ExactInteger] | s.Tuple
)
type ExactScalar = s.Expr | int
type Interval = Literal["0<z<1"]
type Boundary = Literal["upper", "lower"]
type RealRegion = Literal["0<x<1", "-1<x<0", "x>1", "x<-1"]
