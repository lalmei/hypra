"""Direct comparisons on the negative real ray with the pinned HypExp package."""

import pytest
import sympy as s

from hypra import NativeConfig, expand_epsilon

pytestmark = [pytest.mark.real_backend, pytest.mark.usefixtures("require_real_backend")]
eps = s.Symbol("eps")


def test_negative_arguments_and_binomial_reductions_live():
    # Independent symbols distinguish each comparison in one Wolfram session.
    cases = (
        s.hyper([1, eps], [1 + eps], -s.S.Half),
        s.hyper([1, eps], [1 + eps], -2),
        s.hyper([s.S.Half + eps], [], -3),
        s.hyper([s.S.Half + eps, 1 + 2 * eps], [1 + 2 * eps], s.S.Half, evaluate=False),
        s.hyper([1, s.S.Half + eps], [3 * s.S.Half + eps], -1),
    )
    labels = s.symbols("w0:5")
    batch = sum(label * f for label, f in zip(labels, cases, strict=True))
    native = expand_epsilon(batch, eps, 1, config=NativeConfig(use_cache=False))
    upstream = expand_epsilon(batch, eps, 1, backend="wolfram")
    for label in labels:
        difference = native.polynomial.coeff(label) - upstream.polynomial.coeff(label)
        for k in (0, 1):
            assert abs(difference.expand().coeff(eps, k).evalf(40)) < s.Rational(1, 10**32)
