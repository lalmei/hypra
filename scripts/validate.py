#!/usr/bin/env python3
"""Run actual native calculations plus independent analytic validation."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import sympy as s
from validate_integral import derive_j4, derive_moments, numerical_check

from hypra import NativeConfig, backend_info, expand_epsilon


def native_validation():
    eps, w = s.symbols("eps w")
    domain = s.And(w > 0, w < 1)
    p = s.hyper([1, -eps], [1 - 2 * eps], w)
    q = s.hyper([1, -2 * eps], [1 - 2 * eps], 1 - w)
    config = NativeConfig(use_cache=False)
    pr = expand_epsilon(p, eps, 2, domain, config=config)
    qr = expand_epsilon(q, eps, 2, domain, config=config)
    integrand = w ** (-1 - 2 * eps) * (1 - w) ** (-eps) * (p - 1) * (q - 1)
    ir = expand_epsilon(integrand, eps, 3, domain, config=config)
    # The following integration is specific to J4, justified in mathematics.md.
    # Unknown monomials fail; this is not an automatic integral/limit interchange.
    l, m, u, v = s.symbols("L M U V")
    moments = derive_moments()
    substitutions = {s.log(w): l, s.log(1 - w): m, s.polylog(2, w): u, s.polylog(2, 1 - w): v}
    lookup = {
        l * m: moments["log(w)*log(1-w)"],
        l * u: moments["log(w)*Li2(w)"],
        l * m**2: moments["log(w)*log(1-w)^2"],
        m * v: moments["log(1-w)*Li2(1-w)"],
        l**2 * m: moments["log(w)^2*log(1-w)"],
    }
    values = {}
    for k in (2, 3):
        poly = s.expand((ir.coefficient(k) * w).xreplace(substitutions))
        total = s.S.Zero
        for term in s.Add.make_args(poly):
            coefficient, monomial = term.as_coeff_Mul()
            if monomial not in lookup:
                raise ValueError(f"No proved integral for {monomial}")
            total += coefficient * lookup[monomial]
        values[k] = s.simplify(total)
    independent = derive_j4()
    assert s.simplify(values[2] - independent["coefficient_2"]) == 0
    assert s.simplify(values[3] - independent["coefficient_3"]) == 0
    return {
        "P": pr.to_dict(),
        "Q": qr.to_dict(),
        "J4_integrand": ir.to_dict(),
        "J4": values[2] * eps**2 + values[3] * eps**3 + s.Order(eps**4),
        "proven_integrals_checked": 5,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--external", action="store_true", help="Also require experimental Wolfram comparison tests"
    )
    parser.add_argument("--numeric", action="store_true", help="Add secondary quadrature")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    native = native_validation()
    record = {
        "native_backend_validation": "PASSED: actual native Python P/Q/integrand calculations",
        "backend": backend_info(),
        "P": native["P"],
        "Q": native["Q"],
        "J4_integrand": native["J4_integrand"],
        "J4": str(native["J4"]),
        "integration_scope": "J4 only, using five independently derived moments and the documented endpoint proof",
        "external_validation": "not required; not requested",
    }
    if args.numeric:
        record["secondary_quadrature"] = numerical_check()
    code = 0
    if args.external:
        info = backend_info(backend="wolfram")
        record["external_backend"] = info
        if not info["available"]:
            record["external_validation"] = "BLOCKED: " + str(info["error"])
            code = 2
        else:
            import os

            env = dict(os.environ, HYPRA_RUN_REAL="1")
            run = subprocess.run(
                [sys.executable, "-m", "pytest", "-q", "tests/test_real_backend.py"],
                cwd=Path(__file__).resolve().parents[1],
                env=env,
                check=False,
            )
            code = run.returncode
            record["external_validation"] = "passed" if code == 0 else "failed"
    text = json.dumps(record, indent=2) + "\n"
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "validation.json").write_text(text)
        latex = "% Actual native Python output; J4 uses separately proved integration.\n"
        for name in ("P", "Q"):
            latex += "\\[" + name + "=" + native[name]["latex"] + "\\]\n"
        latex += "\\[J_4=" + s.latex(native["J4"]) + "\\]\n"
        (args.output_dir / "validation.tex").write_text(latex)
    print(text, end="")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
