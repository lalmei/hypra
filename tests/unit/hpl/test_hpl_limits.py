import pytest
import sympy as s

from hypra.errors import UnsupportedExpressionError
from hypra.hpl.limits import hpl_limit_at_one
from hypra.symbolic import HPL

z = s.Symbol("z", positive=True)


def test_combined_rational_log_limit_not_separate_divergences():
    assert hpl_limit_at_one((s.polylog(2, z) - s.zeta(2)) / (1 - z) - s.log(1 - z), z) == -1
    assert hpl_limit_at_one(s.log(z) / (1 - z), z) == -1
    assert hpl_limit_at_one(s.polylog(2, z) / z, z) == s.zeta(2)


def test_divergent_and_unsupported_reflection_reject():
    for expr in (s.log(1 - z), HPL((1, -1), z), HPL((0,) * 8 + (1,), z), s.polylog(2, s.sqrt(z))):
        with pytest.raises(UnsupportedExpressionError):
            hpl_limit_at_one(expr, z)
