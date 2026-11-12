"""Optional killable native computation process; exact JSON, no source evaluation."""
from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import asdict

from .backend import BackendProtocolError, BackendTimeoutError, _strict_json_loads
from .errors import IncompleteExpansionError, UnsupportedExpressionError


def run_bounded(request, config):
    settings = asdict(config)
    timeout = settings.pop("timeout_seconds")
    payload = json.dumps({"request": request, "config": settings})
    try:
        completed = subprocess.run(
            [sys.executable, "-m", "hypra.native_worker"], input=payload,
            capture_output=True, text=True, timeout=timeout, check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise BackendTimeoutError(f"Native expansion exceeded {timeout:g} seconds") from exc
    if completed.returncode != 0:
        raise BackendProtocolError("Native computation process failed without a valid result")
    try:
        envelope = _strict_json_loads(completed.stdout)
        if not isinstance(envelope, dict):
            raise BackendProtocolError("Native worker response is not an object")
        if "error" in envelope:
            kind = envelope["error"]
            message = envelope.get("message", "Native worker failed")
            if not isinstance(message, str):
                raise BackendProtocolError("Invalid native worker error message")
            classes = {"UnsupportedExpressionError": UnsupportedExpressionError,
                       "IncompleteExpansionError": IncompleteExpansionError}
            raise classes.get(kind, BackendProtocolError)(message)
        response = envelope["response"]
        if not isinstance(response, dict):
            raise BackendProtocolError("Invalid native worker response")
        return response
    except (KeyError, json.JSONDecodeError) as exc:
        raise BackendProtocolError("Invalid native worker response") from exc


def main():
    from .native import NativeBackend, NativeConfig
    try:
        payload = _strict_json_loads(sys.stdin.read())
        settings = payload["config"]
        settings["timeout_seconds"] = None
        response = NativeBackend(NativeConfig(**settings)).expand(payload["request"])
        result = {"response": response}
    except Exception as exc:  # noqa: BLE001 - child-process exception boundary
        result = {"error": type(exc).__name__, "message": str(exc)}
    print(json.dumps(result))


if __name__ == "__main__":
    main()
