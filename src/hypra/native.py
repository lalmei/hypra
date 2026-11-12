"""Native exact SymPy backend for proved hypergeometric families.

The principal 2F1 family uses an all-order differential-equation/HPL recursion,
not an example lookup table. This is an initial native port, not HypExp parity.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

import sympy as s

from .errors import IncompleteExpansionError, UnsupportedExpressionError
from .planning import _positive, plan_expansion
from .serialization import decode_expr, encode_expr
from .symbolic import HPL, reduce_known_hpl


@dataclass(frozen=True)
class NativeConfig:
    max_internal_order: int = 12
    use_cache: bool = True
    cache_dir: str | None = None
    timeout_seconds: float | None = None
    argument_domain: str = "real-interval"

    def __post_init__(self):
        if self.argument_domain not in ("real-interval", "principal-disk"):
            raise ValueError("argument_domain must be real-interval or principal-disk")
        if self.timeout_seconds is not None and (
            isinstance(self.timeout_seconds, bool)
            or not isinstance(self.timeout_seconds, (int, float))
            or not 0 < self.timeout_seconds < float("inf")
        ):
            raise ValueError("timeout_seconds must be a finite positive number or None")
        if type(self.max_internal_order) is not int or self.max_internal_order < 0:
            raise ValueError("max_internal_order must be a nonnegative integer")


@dataclass(frozen=True)
class NativeInfo:
    available: bool = True
    backend: str = "native-sympy"
    algorithm_version: str = "native-contiguous-nested-sums-v3"
    sympy_version: str = s.__version__
    error: str | None = None

    def as_dict(self):
        return asdict(self)


def _inside_interval(z, assumptions):
    return _positive(z, assumptions) and _positive(1 - z, assumptions)


def _linear(expr, eps, base):
    """Return rational slope if expr=base+slope*eps, else None."""
    slope = s.expand(expr).coeff(eps)
    if slope.is_Rational and s.expand(expr - base - slope * eps) == 0:
        return slope
    return None


def _zero_base_twof1(b, c, z, eps, order):
    # f1=b H_1; fk= -c H_(0,word)+(b-c) H_(1,word).
    result = s.S.One
    words = {(1,): b}
    for k in range(1, order + 1):
        result += eps**k * sum(coef * HPL(word, z) for word, coef in words.items())
        following = {}
        for word, coef in words.items():
            for letter, factor in ((0, -c), (1, b - c)):
                value = coef * factor
                if value:
                    following[(letter, *word)] = value
        words = following
    return reduce_known_hpl(result)


def _li(k, z):
    return -s.log(1 - z) if k == 1 else s.polylog(k, z)


def _half_telescoping(b, z, eps, order):
    t = s.sqrt(z)

    def odd_sum(k):
        return s.Integer(2) ** (k - 1) * (_li(k, t) - _li(k, -t)) / t

    result = (s.log(1 + t) - s.log(1 - t)) / (2 * t)
    for k in range(1, order + 1):
        result += (-1) ** (k - 1) * b**k * (odd_sum(k) - odd_sum(k + 1) / 2) * eps**k
    return result


def _half_even(a, z, eps, order):
    t = s.sqrt(z)
    return sum(
        ((-2 * a) ** k / (2 * s.factorial(k))) * (s.log(1 - t) ** k + s.log(1 + t) ** k) * eps**k
        for k in range(order + 1)
    )


def expand_atom(h, eps, order, assumptions=s.true, argument_domain="real-interval"):
    """Exact coefficient recursion in the explicit real domain 0<z<1."""
    if h.argument == 0:
        return s.S.One
    z = h.argument
    if argument_domain == "principal-disk":
        # This increment admits exact numerical complex points only. Symbolic
        # disk constraints require a separate assumption proof, not a guess.
        if z.free_symbols or not _positive(1-s.Abs(z), s.true) or not (
            z.is_real is False or _positive(z, s.true)
        ):
            raise UnsupportedExpressionError(
                "Principal-disk expansion requires an exact numeric argument with "
                "0 < |z| < 1 off the negative real axis"
            )
    elif not _inside_interval(z, assumptions):
        raise UnsupportedExpressionError(
            "Native hypergeometric expansion requires 0 < argument < 1"
        )
    if len(h.ap) == 1 and not h.bq:
        return s.series((1 - z) ** (-h.ap[0]), eps, 0, order + 1).removeO()
    # Preserve compact classical forms for two common half-integer subfamilies.
    if len(h.ap) == 2 and len(h.bq) == 1:
        upper, lower = tuple(h.ap), h.bq[0]
        if s.S.One in upper:
            other = upper[1] if upper[0] == 1 else upper[0]
            b = _linear(other, eps, s.Rational(1, 2))
            c = _linear(lower, eps, s.Rational(3, 2))
            if b is not None and b == c:
                return _half_telescoping(b, z, eps, order)
        if lower == s.Rational(1, 2):
            for first, second in (upper, upper[::-1]):
                a = _linear(first, eps, s.S.Zero)
                b = _linear(second, eps, s.Rational(1, 2))
                if a is not None and a == b:
                    return _half_even(a, z, eps, order)
    from .canonical import build_canonical, canonical_bases, use_contiguous
    from .contiguous import expand_contiguous
    from .matched_half import half_parameters

    parsed = half_parameters(h, eps)
    if parsed is not None and use_contiguous(*parsed):
        ap, bq = parsed
        ca, cb = canonical_bases(ap, bq)
        coordinate = s.Dummy("argument", positive=True)
        result = expand_contiguous(ap, bq, ca, cb, coordinate, eps, order, build_canonical)
        return reduce_known_hpl(result).xreplace({coordinate: z})
    from .integer import expand_integer_hyper

    integer_result = expand_integer_hyper(h, eps, order)
    if integer_result is not None:
        return integer_result
    from .matched_half import expand_matched_half

    matched_result = expand_matched_half(h, eps, order)
    if matched_result is not None:
        return matched_result
    if len(h.ap) != 2 or len(h.bq) != 1:
        from .half_system import canonical_half_series
        from .matched_half import half_parameters

        parsed = half_parameters(h, eps)
        if parsed is not None:
            ap, bq = parsed
            result = canonical_half_series(
                [a for a, _ in ap], [a for _, a in ap],
                [b for b, _ in bq], [b for _, b in bq], z, eps, order,
            )
            if result is not None:
                return result
        raise UnsupportedExpressionError(
            "Unsupported noninteger parameter pattern or shifts"
        )
    upper, lower = tuple(h.ap), h.bq[0]
    if s.S.One in upper:
        other = upper[1] if upper[0] == 1 else upper[0]
        b, c = _linear(other, eps, s.S.Zero), _linear(lower, eps, s.S.One)
        if b is not None and c is not None:
            return _zero_base_twof1(b, c, z, eps, order)
        b = _linear(other, eps, s.Rational(1, 2))
        c = _linear(lower, eps, s.Rational(3, 2))
        if b is not None and b == c:
            return _half_telescoping(b, z, eps, order)
    if lower == s.Rational(1, 2):
        for first, second in (upper, upper[::-1]):
            a = _linear(first, eps, s.S.Zero)
            b = _linear(second, eps, s.Rational(1, 2))
            if a is not None and a == b:
                return _half_even(a, z, eps, order)
    from .twof1_shift import expand_shifted_twof1

    shifted = expand_shifted_twof1(h, eps, order)
    if shifted is not None:
        return shifted
    from .half_system import canonical_half_series
    from .matched_half import half_parameters

    parsed = half_parameters(h, eps)
    if parsed is not None:
        ap, bq = parsed
        canonical = canonical_half_series(
            [a for a, _ in ap],
            [a for _, a in ap],
            [b for b, _ in bq],
            [b for _, b in bq],
            z,
            eps,
            order,
        )
        if canonical is not None:
            return reduce_known_hpl(canonical)
    raise UnsupportedExpressionError(
        "Native port does not yet support this parameter family; no external fallback was used"
    )


class NativeBackend:
    def __init__(self, config=None):
        self.config = config or NativeConfig()
        if not isinstance(self.config, NativeConfig):
            raise TypeError("Native backend requires NativeConfig")

    def info(self):
        return NativeInfo()

    def _fingerprint(self):
        digest = hashlib.sha256()
        for name in (
            "native.py",
            "planning.py",
            "serialization.py",
            "symbolic.py",
            "integer.py",
            "nested_sums.py",
            "hpl_algebra.py",
            "colored_sums.py",
            "matched_half.py",
            "half_system.py",
            "twof1_shift.py",
            "endpoints.py",
            "canonical.py",
            "contiguous.py",
        ):
            digest.update(Path(__file__).with_name(name).read_bytes())
        return digest.hexdigest()

    def expand(self, request):
        if self.config.timeout_seconds is not None:
            from .native_worker import run_bounded

            return run_bounded(request, self.config)
        return self._expand(request)

    def _expand(self, request):
        expr, eps = decode_expr(request["expression"]), decode_expr(request["epsilon"])
        assumptions = decode_expr(request["assumptions"])
        plan = plan_expansion(expr, eps, request["through"], assumptions)
        order = request["internal_order"]
        if type(order) is not int or order < plan.internal_order:
            raise IncompleteExpansionError("Insufficient native internal order")
        from .canonical import canonical_bases, use_contiguous
        from .contiguous import plan_contiguous
        from .integer import expansion_depth
        from .matched_half import half_parameters

        def atom_depth(h):
            parsed = half_parameters(h, eps)
            if parsed is not None and use_contiguous(*parsed):
                return order + plan_contiguous(*parsed, *canonical_bases(*parsed)).extra_order
            return expansion_depth(h, eps, order)

        canonical_order = max(
            [order, *(atom_depth(h) for h in expr.atoms(s.hyper))]
        )
        if canonical_order > self.config.max_internal_order:
            raise UnsupportedExpressionError(
                f"Internal order {canonical_order} exceeds configured native limit {self.config.max_internal_order}"
            )
        meta = self.info().as_dict() | {
            "source_fingerprint": self._fingerprint(),
            "canonical_order_max": canonical_order,
            "argument_domain": self.config.argument_domain,
        }
        payload = json.dumps(
            {"request": request, "backend": meta}, sort_keys=True, separators=(",", ":")
        )
        key = hashlib.sha256(payload.encode()).hexdigest()
        root = (
            Path(self.config.cache_dir).expanduser()
            if self.config.cache_dir
            else Path.home() / ".cache" / "hypra" / "native"
        )
        cache = root / (key + ".json")
        if self.config.use_cache and cache.exists():
            try:
                from .backend import _strict_json_loads

                envelope = _strict_json_loads(cache.read_text())
                encoded = json.dumps(envelope["response"], sort_keys=True, separators=(",", ":"))
                if (
                    envelope["key"] == key
                    and envelope["checksum"] == hashlib.sha256(encoded.encode()).hexdigest()
                ):
                    response = envelope["response"]
                    response["cache_hit"] = True
                    return response
            except (OSError, ValueError, KeyError, TypeError):
                pass
        replacements = {h: expand_atom(h, eps, order, assumptions, self.config.argument_domain)
                        for h in expr.atoms(s.hyper)}
        expanded = expr.xreplace(replacements)
        try:
            if request["through"] < plan.valuation_lower_bound:
                polynomial = s.S.Zero
            else:
                shifted = expanded * eps**plan.pole_bound
                series = s.series(shifted, eps, 0, request["through"] + plan.pole_bound + 1)
                polynomial = s.expand(series.removeO() / eps**plan.pole_bound)
        except (NotImplementedError, ValueError, TypeError, RecursionError) as exc:
            raise IncompleteExpansionError(
                "SymPy could not complete the final Laurent expansion"
            ) from exc
        if polynomial.has(s.hyper, s.Derivative, s.Subs, s.nan, s.zoo, s.oo, -s.oo):
            raise IncompleteExpansionError(
                "Native expansion contains unresolved or singular coefficients"
            )
        response = {
            "result": encode_expr(polynomial),
            "complete": True,
            "achieved_through": request["through"],
            "remainder_order": request["through"] + 1,
            "backend": meta,
            "provenance": meta,
            "order_basis": "proved-native-recurrence-and-conservative-pole-budget",
            "cache_hit": False,
        }
        if self.config.use_cache:
            temporary = None
            try:
                root.mkdir(parents=True, exist_ok=True)
                encoded = json.dumps(response, sort_keys=True, separators=(",", ":"))
                with tempfile.NamedTemporaryFile("w", dir=root, delete=False) as stream:
                    temporary = Path(stream.name)
                    json.dump(
                        {
                            "key": key,
                            "checksum": hashlib.sha256(encoded.encode()).hexdigest(),
                            "response": response,
                        },
                        stream,
                    )
                os.replace(temporary, cache)
            except OSError:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
        return response
