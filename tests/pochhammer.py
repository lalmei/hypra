"""Exact reference coefficients for the nested-sum tests."""

import sympy as s


def pochhammer_taylor(upper, lower, regulator, n, order):
    """Expand the defining coefficient of z**n through the given regulator order."""
    numerator = s.Poly(
        s.prod(s.rf(base + slope * regulator, n) for base, slope in upper), regulator
    )
    denominator = s.Poly(
        s.factorial(n) * s.prod(s.rf(base + slope * regulator, n) for base, slope in lower),
        regulator,
    )
    constant = denominator.nth(0)
    assert constant != 0

    # Divide the finite Pochhammer products as power series at regulator zero.
    coefficients = []
    for degree in range(order + 1):
        previous = sum(denominator.nth(k) * coefficients[degree - k] for k in range(1, degree + 1))
        coefficients.append((numerator.nth(degree) - previous) / constant)
    return sum(value * regulator**degree for degree, value in enumerate(coefficients))
