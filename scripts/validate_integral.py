#!/usr/bin/env python3
"""Independent analytic validation of the P, Q and J4 acceptance example.

This derives the necessary logarithmic moments, rather than importing the
NNLO report's lookup table. It is not an implementation of the Wolfram bridge
and does not establish that a real HypExp installation works.
"""

from __future__ import annotations

import argparse
import json

import sympy as s

EPS = s.Symbol("epsilon")
W = s.Symbol("w", positive=True)


def derive_moments() -> dict[str, s.Expr]:
    """Derive five moments by convergent sums, beta derivatives and reflection.

    Each moment is integrated with measure dw/w over (0,1). See the mathematical
    documentation for the justified exchanges and the beta coefficient proof.
    """
    n = s.Symbol("n", integer=True, positive=True)
    a, t = s.symbols("a t")
    log_gamma_ratio = s.expand(
        sum((-1) ** k * s.zeta(k) * (a**k + t**k - (a + t) ** k) / k for k in range(2, 5))
    )
    # Beta(a, 1+t) = exp(log_gamma_ratio)/a. Taking d_a d_t^2
    # means twice the a^2*t^2 coefficient before division by a.
    beta_mixed = s.expand(log_gamma_ratio + log_gamma_ratio**2 / 2)
    b = s.simplify(2 * beta_mixed.coeff(a, 2).coeff(t, 2))
    zeta2_from_series = s.summation(n**-2, (n, 1, s.oo))
    # Li2(1-w) = zeta(2)-Li2(w)-log(w)log(1-w), and
    # integral Li2(w)log(1-w) dw/w = -Li2(1)^2/2.
    c = s.simplify(-(zeta2_from_series**2) / 2 - b)
    return {
        "log(w)*log(1-w)": s.summation(n**-3, (n, 1, s.oo)),
        "log(w)*Li2(w)": -s.summation(n**-4, (n, 1, s.oo)),
        "log(w)*log(1-w)^2": b,
        "log(1-w)*Li2(1-w)": c,
        "log(w)^2*log(1-w)": -2 * s.summation(n**-4, (n, 1, s.oo)),
    }


def derive_j4() -> dict[str, s.Expr | dict[str, s.Expr]]:
    """Assemble the Taylor coefficients from derived, exact moments."""
    l, m, u, v = s.symbols("L M U V")
    p = EPS * m + EPS**2 * (-2 * u - m**2 / 2)
    q = 2 * EPS * l - 4 * EPS**2 * v
    weight = 1 + EPS * (-2 * l - m)
    integrand = s.expand(p * q * weight)
    moments = derive_moments()
    monomial_integrals = {
        l * m: moments["log(w)*log(1-w)"],
        l * u: moments["log(w)*Li2(w)"],
        l * m**2: moments["log(w)*log(1-w)^2"],
        m * v: moments["log(1-w)*Li2(1-w)"],
        l**2 * m: moments["log(w)^2*log(1-w)"],
    }

    def integrate_polynomial(poly: s.Expr) -> s.Expr:
        result = s.S.Zero
        for term in s.Add.make_args(s.expand(poly)):
            coefficient, monomial = term.as_coeff_Mul()
            if monomial not in monomial_integrals:
                raise ValueError(f"Unproved logarithmic moment: {monomial}")
            result += coefficient * monomial_integrals[monomial]
        return s.simplify(result)

    c2 = integrate_polynomial(integrand.coeff(EPS, 2))
    c3 = integrate_polynomial(integrand.coeff(EPS, 3))
    substitutions = {
        l: s.log(W),
        m: s.log(1 - W),
        u: s.polylog(2, W),
        v: s.polylog(2, 1 - W),
    }
    return {
        "P": p.subs(substitutions) + 1 + s.Order(EPS**3, EPS),
        "Q": q.subs(substitutions) + 1 + s.Order(EPS**3, EPS),
        "J4": c2 * EPS**2 + c3 * EPS**3 + s.Order(EPS**4, EPS),
        "coefficient_2": c2,
        "coefficient_3": c3,
        "integrand_coefficient_2": integrand.coeff(EPS, 2).subs(substitutions) / W,
        "integrand_coefficient_3": integrand.coeff(EPS, 3).subs(substitutions) / W,
        "moments": moments,
    }


def numerical_check(dps: int = 40) -> list[dict[str, str]]:
    """Secondary quadrature, with no endpoint cutoffs or fitted coefficients.

    Negative epsilon makes endpoint limits zero. Contiguous hypergeometric
    representations avoid subtraction of nearly equal P and 1 (or Q and 1).
    The residual is expected to be O(epsilon**4), not exactly zero.
    """
    import mpmath as mp

    rows = []
    with mp.workdps(dps):
        for eps in map(mp.mpf, ("-0.01", "-0.005", "-0.0025")):

            def integrand(w, eps=eps):
                if w == 0 or w == 1:
                    return mp.mpf(0)
                pm1 = -eps * w / (1 - 2 * eps) * mp.hyp2f1(1, 1 - eps, 2 - 2 * eps, w)
                qm1 = (
                    -2
                    * eps
                    * (1 - w)
                    / (1 - 2 * eps)
                    * mp.hyp2f1(1, 1 - 2 * eps, 2 - 2 * eps, 1 - w)
                )
                return w ** (-1 - 2 * eps) * (1 - w) ** (-eps) * pm1 * qm1

            value = mp.quad(integrand, [0, mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("0.9"), 1])
            derived = derive_j4()
            c2 = mp.mpf(str(s.N(derived["coefficient_2"], dps)))
            c3 = mp.mpf(str(s.N(derived["coefficient_3"], dps)))
            residual = value - c2 * eps**2 - c3 * eps**3
            rows.append(
                {
                    "epsilon": str(eps),
                    "J4_quadrature": mp.nstr(value, dps - 5),
                    "residual_over_epsilon4": mp.nstr(residual / eps**4, 15),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numeric", action="store_true", help="Add secondary mpmath quadrature")
    parser.add_argument("--latex", action="store_true", help="Print derivation outputs as LaTeX")
    args = parser.parse_args()
    derived = derive_j4()
    if args.latex:
        for name in ("P", "Q", "J4"):
            print(r"\[" + name + " = " + s.latex(derived[name]) + r"\]")
    else:
        result = {
            "kind": "independent analytic derivation; not a Wolfram backend validation",
            "domain": "0 < w < 1, epsilon near 0; real principal branches",
            "P": str(derived["P"]),
            "Q": str(derived["Q"]),
            "J4": str(derived["J4"]),
            "derived_moments": {key: str(value) for key, value in derived["moments"].items()},
        }
        if args.numeric:
            result["secondary_quadrature"] = numerical_check()
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
