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
from typing import Any, Literal

import sympy as s
from sympy.logic.boolalg import Boolean

from .._json import strict_json_loads
from ..errors import IncompleteExpansionError, UnsupportedExpressionError
from ..serialization import decode_expr, encode_expr
from .atom import expand_atom
from .planning import plan_expansion


@dataclass(frozen=True)
class NativeConfig:
    """Settings for the native SymPy expansion engine.

    Attributes:
        max_internal_order: Maximum internal power after accounting for poles
            and reductions.
        use_cache: Reuse exact expansions from the local disk cache.
        cache_dir: Override the default native cache directory.
        timeout_seconds: Optional worker-process timeout in seconds.
        expansion_library: Optional SQLite file of registered symbolic families.
        argument_domain: Domain contract for interior hypergeometric arguments.
        unit_endpoint: Require convergence at z=1 or allow supported
            meromorphic parameter continuation.
        unit_excess_shift: One-based index of a lower parameter to shift at a
            fixed excess pole in meromorphic mode.
        unit_shift_singular_lowers: Add epsilon to each fixed nonpositive integer
            lower parameter at z=1 before cancellation in meromorphic mode.
    """

    max_internal_order: int = 12
    use_cache: bool = True
    cache_dir: str | None = None
    timeout_seconds: float | None = None
    argument_domain: Literal["real-interval", "principal-disk", "principal-cut-plane"] = (
        "real-interval"
    )
    unit_endpoint: Literal["convergent", "meromorphic"] = "convergent"
    unit_excess_shift: int | None = None
    unit_shift_singular_lowers: bool = False
    expansion_library: str | None = None

    def __post_init__(self) -> None:
        if self.expansion_library is not None and not isinstance(self.expansion_library, str):
            raise TypeError("expansion_library must be a path string")
        if self.unit_endpoint not in ("convergent", "meromorphic"):
            raise ValueError("unit_endpoint must be convergent or meromorphic")
        if type(self.unit_shift_singular_lowers) is not bool or (
            self.unit_shift_singular_lowers and self.unit_endpoint != "meromorphic"
        ):
            raise ValueError("unit_shift_singular_lowers requires a Boolean and meromorphic mode")
        if self.unit_excess_shift is not None and (
            self.unit_endpoint != "meromorphic"
            or type(self.unit_excess_shift) is not int
            or self.unit_excess_shift < 1
        ):
            raise ValueError("unit_excess_shift requires meromorphic mode and a positive index")
        if self.argument_domain not in ("real-interval", "principal-disk", "principal-cut-plane"):
            raise ValueError(
                "argument_domain must be real-interval, principal-disk or principal-cut-plane"
            )
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
    """Availability and version details for the native engine."""

    available: bool = True
    backend: str = "native-sympy"
    algorithm_version: str = "native-contiguous-nested-sums-v3"
    sympy_version: str = s.__version__
    error: str | None = None

    def as_dict(self) -> dict[str, Any]:
        """Return availability and version details as a JSON-compatible mapping."""
        return asdict(self)


class NativeBackend:
    """Execute exact AST expansion requests with the native recurrence engine.

    Args:
        config: Native settings. Omit to use the defaults.
    """

    def __init__(self, config: NativeConfig | None = None) -> None:
        self.config = config or NativeConfig()
        if not isinstance(self.config, NativeConfig):
            raise TypeError("Native backend requires NativeConfig")

    def info(self) -> NativeInfo:
        """Return availability and version information."""
        return NativeInfo()

    def _fingerprint(self) -> str:
        digest = hashlib.sha256()
        package = Path(__file__).resolve().parents[1]
        for name in (
            "native/backend.py",
            "expansion_library.py",
            "native/atom.py",
            "native/planning.py",
            "serialization.py",
            "symbolic.py",
            "native/integer.py",
            "native/nested_sums.py",
            "hpl/algebra.py",
            "hpl/limits.py",
            "hpl/transforms.py",
            "hpl/endpoints.py",
            "native/colored_sums.py",
            "native/matched_half.py",
            "native/half_system.py",
            "native/twof1_shift.py",
            "native/endpoints.py",
            "native/unit_limit.py",
            "hpl/series_one.py",
            "hpl/inversion.py",
            "hpl/mobius.py",
            "native/canonical.py",
            "native/contiguous.py",
        ):
            digest.update((package / name).read_bytes())
        return digest.hexdigest()

    def expand(self, request: dict[str, Any]) -> dict[str, Any]:
        """Execute a protocol expansion request and return the exact result AST.

        Args:
            request: Mapping containing encoded input and the requested orders.

        Returns:
            Protocol response with the result, achieved order, and provenance.
        """
        if self.config.timeout_seconds is not None:
            from .worker import run_bounded

            return run_bounded(request, self.config)
        return self._expand(request)

    def _expand(self, request: dict[str, Any]) -> dict[str, Any]:
        expr, eps = decode_expr(request["expression"]), decode_expr(request["epsilon"])
        assumptions = decode_expr(request["assumptions"])
        if not isinstance(expr, s.Expr) or not isinstance(eps, s.Symbol):
            raise TypeError("Expected a SymPy expression and regulator Symbol")
        if not isinstance(assumptions, Boolean):
            raise TypeError("Assumptions must be a SymPy Boolean")
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

        canonical_order = max([order, *(atom_depth(h) for h in expr.atoms(s.hyper))])
        if canonical_order > self.config.max_internal_order:
            raise UnsupportedExpressionError(
                f"Internal order {canonical_order} exceeds configured native limit {self.config.max_internal_order}"
            )
        meta = self.info().as_dict() | {
            "source_fingerprint": self._fingerprint(),
            "canonical_order_max": canonical_order,
            "argument_domain": self.config.argument_domain,
            "unit_endpoint": self.config.unit_endpoint,
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
                envelope = strict_json_loads(cache.read_text())
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

        library = None
        if (
            self.config.expansion_library is not None
            and self.config.argument_domain == "real-interval"
        ):
            from ..expansion_library import ExpansionLibrary

            library = ExpansionLibrary(self.config.expansion_library)

        def checked_atom(h):
            from .atom import _inside_interval

            result = None
            if library is not None and _inside_interval(h.argument, assumptions):
                result = library._lookup_atom(h, eps, order)
            if result is None:
                result = expand_atom(h, eps, order, assumptions, self.config.argument_domain)
            slope_symbols = (
                set().union(*(s.diff(parameter, eps).free_symbols for parameter in (*h.ap, *h.bq)))
                if h.ap or h.bq
                else set()
            )
            if slope_symbols:
                # Each nonpositive integer lower base contributes at most one
                # genuine inverse slope factor. Other parameter denominators
                # from contiguous operators must cancel exactly.
                source_poles = s.prod(
                    s.diff(parameter, eps)
                    for parameter in h.bq
                    if parameter.subs(eps, 0).is_Integer and parameter.subs(eps, 0) <= 0
                )
                from .contiguous import _hpl_pieces, _normalize_coefficient

                result = _normalize_coefficient(result)
                if any(
                    s.cancel(source_poles * coefficient).is_polynomial(*slope_symbols) is not True
                    for coefficient in _hpl_pieces(result).values()
                ):
                    raise IncompleteExpansionError(
                        "Symbolic regulator slopes left an unresolved parameter denominator"
                    )
                result = s.expand(result)
            return result

        replacements = {h: checked_atom(h) for h in expr.atoms(s.hyper)}
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
