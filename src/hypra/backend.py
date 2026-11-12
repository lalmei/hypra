"""Fail-closed JSON/file subprocess transport for external Wolfram packages.

No shell, eval, arbitrary expression strings, or Mathematica source generated
from user expressions is used. The bundled bridge interprets a finite AST.
"""

from __future__ import annotations

import hashlib
import json
import os
import signal
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .config import BackendConfig

PROTOCOL_VERSION = 1
MAX_RESPONSE_BYTES = 16 * 1024 * 1024


class BackendError(RuntimeError):
    """A backend failure; no partial mathematical result is accepted."""


class BackendUnavailableError(BackendError):
    """The external licensed backend or one of its packages is unavailable."""


class BackendTimeoutError(BackendError):
    """The backend exceeded its configured wall-clock limit."""


class BackendProtocolError(BackendError):
    """The backend did not return a complete, valid protocol response."""


@dataclass(frozen=True)
class BackendInfo:
    available: bool
    executable: str | None
    wolfram_version: str | None = None
    system_id: str | None = None
    hypexp_version: str | None = None
    hpl_version: str | None = None
    hypexp_path: str | None = None
    hpl_path: str | None = None
    fingerprint: str | None = None
    hypexp_sha256: str | None = None
    hpl_sha256: str | None = None
    error: str | None = None
    bridge_version: int = PROTOCOL_VERSION

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _strict_json_loads(text: str) -> Any:
    """Reject duplicate object keys and JSON's nonstandard nonfinite constants."""

    def unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> Any:
        raise ValueError(f"Nonfinite JSON constant: {value}")

    return json.loads(text, object_pairs_hook=unique_pairs, parse_constant=reject_constant)


def _file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _package_fingerprint(path: str | None) -> str | None:
    """Hash package source siblings too, so helper-file changes invalidate cache."""
    if not path:
        return None
    root = Path(path).expanduser()
    if not root.exists():
        return None
    if root.is_file():
        root = root.parent
        if root.name == "Kernel":
            root = root.parent
    digest = hashlib.sha256()
    files = sorted(
        p for p in root.rglob("*") if p.is_file() and p.suffix in {".m", ".wl", ".mx", ".wls"}
    )
    for file in files:
        digest.update(str(file.relative_to(root)).encode())
        digest.update(bytes.fromhex(_file_hash(file)))
    return digest.hexdigest()


class WolframBackend:
    """One-shot kernel execution; cache keys include full requests and versions."""

    def __init__(self, config: BackendConfig | None = None):
        self.config = config or BackendConfig.from_env()
        self.bridge_path = Path(__file__).with_name("wolfram_bridge.wls")

    def _run(self, request: dict[str, Any]) -> dict[str, Any]:
        executable = self.config.resolved_executable()
        if not executable:
            raise BackendUnavailableError(
                "Wolfram executable not found. Install/license Wolfram and configure "
                "HYPRA_WOLFRAM_EXECUTABLE; HypExp and HPL are separate external dependencies."
            )
        if not self.bridge_path.is_file():
            raise BackendUnavailableError("Bundled Wolfram JSON bridge is missing")
        payload = dict(request)
        payload["protocol"] = PROTOCOL_VERSION
        payload["packages"] = {
            "hypexp": self.config.hypexp_path,
            "hpl": self.config.hpl_path,
        }
        with tempfile.TemporaryDirectory(prefix="hypra-") as tmp:
            request_path, response_path = Path(tmp) / "request.json", Path(tmp) / "response.json"
            request_path.write_bytes(_canonical(payload))
            request_path.chmod(0o600)
            mode = "-file" if "wolframscript" in Path(executable).name.lower() else "-script"
            command = [
                executable,
                mode,
                str(self.bridge_path),
                str(request_path),
                str(response_path),
            ]
            try:
                process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    shell=False,
                    start_new_session=(os.name == "posix"),
                )
            except OSError as exc:
                raise BackendUnavailableError(f"Cannot execute Wolfram backend: {exc}") from exc
            try:
                process.communicate(timeout=self.config.timeout)
            except subprocess.TimeoutExpired as exc:
                if os.name == "posix":
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                else:
                    process.kill()
                process.communicate()
                raise BackendTimeoutError(
                    f"Wolfram backend exceeded {self.config.timeout:g} seconds; no result accepted"
                ) from exc
            if process.returncode:
                raise BackendError(
                    f"Wolfram exited with status {process.returncode}; no result accepted"
                )
            if not response_path.is_file() or response_path.stat().st_size > MAX_RESPONSE_BYTES:
                raise BackendProtocolError("Missing or oversized Wolfram response")
            try:
                response = _strict_json_loads(response_path.read_text())
            except (OSError, UnicodeError, ValueError) as exc:
                raise BackendProtocolError("Wolfram response is not valid JSON") from exc
        if (
            not isinstance(response, dict)
            or type(response.get("protocol")) is not int
            or response.get("protocol") != PROTOCOL_VERSION
        ):
            raise BackendProtocolError("Unexpected Wolfram protocol")
        if response.get("status") != "ok":
            detail = response.get("error", "unspecified backend error")
            raise BackendError(f"Wolfram rejected the request: {detail}")
        return response

    def info(self) -> BackendInfo:
        executable = self.config.resolved_executable()
        if executable is None:
            return BackendInfo(False, None, error="Wolfram executable not found")
        try:
            response = self._run({"operation": "info"})
            meta = response.get("backend")
            if not isinstance(meta, dict) or not isinstance(meta.get("wolfram_version"), str):
                raise BackendProtocolError("Backend version metadata missing")
            for key in ("hypexp_path", "hpl_path", "hypexp_sha256", "hpl_sha256"):
                if not isinstance(meta.get(key), str) or not meta[key]:
                    raise BackendProtocolError(f"Backend fingerprint metadata missing: {key}")
            fingerprint_inputs = {
                "reported": meta,
                "bridge_sha256": _file_hash(self.bridge_path),
                "executable": executable,
                "executable_sha256": _file_hash(Path(executable)),
                "hypexp_tree": _package_fingerprint(meta["hypexp_path"]),
                "hpl_tree": _package_fingerprint(meta["hpl_path"]),
            }
            return BackendInfo(
                True,
                executable,
                wolfram_version=meta["wolfram_version"],
                system_id=meta.get("system_id"),
                hypexp_version=meta.get("hypexp_version"),
                hpl_version=meta.get("hpl_version"),
                hypexp_path=meta["hypexp_path"],
                hpl_path=meta["hpl_path"],
                fingerprint=hashlib.sha256(_canonical(fingerprint_inputs)).hexdigest(),
                hypexp_sha256=meta["hypexp_sha256"],
                hpl_sha256=meta["hpl_sha256"],
            )
        except (BackendError, OSError, ValueError) as exc:
            return BackendInfo(False, executable, error=str(exc))

    @staticmethod
    def _validate_expansion(response: dict[str, Any], request: dict[str, Any]) -> None:
        through = request["through"]
        if (
            not isinstance(response, dict)
            or type(response.get("protocol")) is not int
            or response.get("protocol") != PROTOCOL_VERSION
            or response.get("status") != "ok"
        ):
            raise BackendProtocolError("Invalid expansion envelope")
        if response.get("complete") is not True:
            raise BackendProtocolError("Backend expansion is incomplete")
        if (
            type(response.get("achieved_through")) is not int
            or type(response.get("remainder_order")) is not int
            or response.get("achieved_through") != through
            or response.get("remainder_order") != through + 1
        ):
            raise BackendProtocolError("Backend truncation order does not match request")
        if not isinstance(response.get("result"), dict):
            raise BackendProtocolError("Backend result AST is missing")
        if response.get("order_basis") != "hypra-pole-bound-plus-hypexp-series":
            raise BackendProtocolError("Backend omitted order provenance")

    def expand(self, request: dict[str, Any]) -> dict[str, Any]:
        """Expand an AST request; caller validates expression family and pole bound.

        ``internal_order`` must cover all series depth required by prefactors.
        The bridge's order annotation records the algorithm's truncation contract,
        not an independent proof of HypExp's correctness or endpoint integrability.
        """
        for key in ("through", "internal_order"):
            if isinstance(request.get(key), bool) or not isinstance(request.get(key), int):
                raise TypeError(f"{key} must be an integer")
        if request["internal_order"] < max(0, request["through"]):
            raise ValueError("internal_order must cover nonnegative requested expansion order")
        info = self.info()
        if not info.available:
            raise BackendUnavailableError(info.error or "External backend unavailable")
        payload = dict(request, operation="expand")
        digest = hashlib.sha256(
            _canonical(
                {
                    "request": payload,
                    "backend": info.as_dict(),
                    "protocol": PROTOCOL_VERSION,
                }
            )
        ).hexdigest()
        cache_path = self.config.resolved_cache_dir() / f"{digest}.json"
        if self.config.use_cache and cache_path.is_file():
            try:
                if cache_path.stat().st_size > MAX_RESPONSE_BYTES:
                    raise ValueError("oversized cache")
                envelope = _strict_json_loads(cache_path.read_text())
                if not isinstance(envelope, dict):
                    raise TypeError("invalid cache envelope")
                if envelope.get("key") != digest:
                    raise ValueError("cache identity mismatch")
                response = envelope["response"]
                if (
                    envelope.get("response_sha256")
                    != hashlib.sha256(_canonical(response)).hexdigest()
                ):
                    raise ValueError("cache checksum mismatch")
                self._validate_expansion(response, payload)
                response["cache_hit"] = True
                response["provenance"] = info.as_dict()
                return response
            except (ValueError, KeyError, TypeError, OSError, BackendError):
                # A corrupt entry is never used as a mathematical result.
                pass
        response = self._run(payload)
        self._validate_expansion(response, payload)
        reported = response.get("backend", {})
        if not isinstance(reported, dict):
            raise BackendProtocolError("Backend identity metadata missing")
        for key in (
            "wolfram_version",
            "system_id",
            "hypexp_path",
            "hpl_path",
            "hypexp_sha256",
            "hpl_sha256",
        ):
            if reported.get(key) != getattr(info, key):
                raise BackendProtocolError(f"Backend identity changed during expansion: {key}")
        response["provenance"] = info.as_dict()
        response["cache_hit"] = False
        if self.config.use_cache:
            # Cache write failure does not invalidate an otherwise verified transport result.
            try:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                with tempfile.NamedTemporaryFile("wb", dir=cache_path.parent, delete=False) as file:
                    temporary = Path(file.name)
                    file.write(
                        _canonical(
                            {
                                "key": digest,
                                "response": response,
                                "response_sha256": hashlib.sha256(_canonical(response)).hexdigest(),
                            }
                        )
                    )
                os.replace(temporary, cache_path)
            except OSError:
                if "temporary" in locals():
                    temporary.unlink(missing_ok=True)
        return response


def backend_info(config: BackendConfig | None = None) -> dict[str, Any]:
    """Return a serializable health/version record, including unavailable states."""
    return WolframBackend(config).info().as_dict()
