"""Offline transport tests. The executable here is a Python stub, NOT Wolfram."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

from hypra.backend import (
    BackendError,
    BackendProtocolError,
    BackendTimeoutError,
    BackendUnavailableError,
    WolframBackend,
    _package_fingerprint,
    backend_info,
)
from hypra.config import BackendConfig


@pytest.fixture
def stub_backend(tmp_path):
    """A protocol stub exercising real subprocess transport only."""
    package_dir = tmp_path / "packages"
    package_dir.mkdir()
    hp, hl = package_dir / "HypExp.m", package_dir / "HPL.m"
    hp.write_text("(* offline stub, no mathematical implementation *)")
    hl.write_text("(* offline stub, no mathematical implementation *)")
    meta = {
        "wolfram_version": "MOCK-NOT-WOLFRAM",
        "system_id": "MOCK",
        "hypexp_version": "MOCK",
        "hpl_version": "MOCK",
        "hypexp_path": str(hp),
        "hpl_path": str(hl),
        "hypexp_sha256": hashlib.sha256(hp.read_bytes()).hexdigest(),
        "hpl_sha256": hashlib.sha256(hl.read_bytes()).hexdigest(),
    }
    executable = tmp_path / "wolfram stub ; harmless"
    executable.write_text(
        f"#!{sys.executable}\n"
        "import json, pathlib, sys, time\n"
        f"meta = {meta!r}\n"
        "request=json.loads(pathlib.Path(sys.argv[-2]).read_text())\n"
        "mode=request.get('stub_mode')\n"
        "output=pathlib.Path(sys.argv[-1])\n"
        "if mode == 'timeout': time.sleep(10)\n"
        "if mode == 'nonzero': sys.exit(7)\n"
        "if mode == 'missing': sys.exit(0)\n"
        "if mode == 'malformed': output.write_text('not json'); sys.exit(0)\n"
        "response={'protocol':1,'status':'ok','backend':meta}\n"
        "if request.get('operation') == 'expand':\n"
        " response.update(result={'head':'Integer','value':'0'},complete=True,"
        "achieved_through=request['through'],remainder_order=request['through']+1,"
        "order_basis='hypra-pole-bound-plus-hypexp-series')\n"
        "if mode == 'incomplete': response['complete']=False\n"
        "if mode == 'wrong_order': response['achieved_through']=-99\n"
        "if mode == 'error': response={'protocol':1,'status':'error','error':'unsupported'}\n"
        "if mode == 'changed_version': response['backend']['wolfram_version']='other'\n"
        "encoded = json.dumps(response)\n"
        "if mode == 'duplicate': encoded=encoded.replace('\"protocol\": 1', '\"protocol\":1,\"protocol\":1')\n"
        "if mode == 'nonfinite': encoded=encoded.replace('\"value\": \"0\"', '\"value\": NaN')\n"
        "if mode == 'boolean_protocol': encoded=encoded.replace('\"protocol\": 1', '\"protocol\": true')\n"
        "output.write_text(encoded)\n"
    )
    executable.chmod(0o700)
    return WolframBackend(
        BackendConfig(
            executable=str(executable),
            hypexp_path=str(hp),
            hpl_path=str(hl),
            timeout=2,
            cache_dir=str(tmp_path / "cache"),
        )
    )


def request(**extra):
    return {
        "expression": {"head": "Integer", "value": "1"},
        "epsilon": {"head": "Symbol", "name": "eps"},
        "assumptions": {"head": "True"},
        "through": 3,
        "internal_order": 3,
        **extra,
    }


def test_missing_backend_is_honest():
    config = BackendConfig(executable="/definitely/not/a/wolfram/executable")
    assert backend_info(config)["available"] is False
    with pytest.raises(BackendUnavailableError):
        WolframBackend(config).expand(request())


def test_stub_health_has_fingerprint(stub_backend):
    info = stub_backend.info()
    assert info.available
    assert info.wolfram_version == "MOCK-NOT-WOLFRAM"
    assert len(info.fingerprint) == 64
    assert len(info.hypexp_sha256) == 64


def test_stub_roundtrip_and_cache(stub_backend):
    result = stub_backend.expand(request())
    assert result["complete"] is True
    assert result["cache_hit"] is False
    assert stub_backend.expand(request())["cache_hit"] is True
    assert stub_backend.expand(request(through=4, internal_order=4))["cache_hit"] is False
    assert stub_backend.expand(request(assumptions={"head": "False"}))["cache_hit"] is False


def test_source_change_invalidates_cache(stub_backend):
    stub_backend.expand(request())
    sibling = Path(stub_backend.config.hypexp_path).with_name("helper.m")
    sibling.write_text("(* different helper *)")
    assert stub_backend.expand(request())["cache_hit"] is False


def test_corrupt_cache_is_not_used(stub_backend):
    stub_backend.expand(request())
    cache = next(stub_backend.config.resolved_cache_dir().glob("*.json"))
    envelope = json.loads(cache.read_text())
    envelope["response"]["result"] = {"head": "Integer", "value": "999"}
    cache.write_text(json.dumps(envelope))
    response = stub_backend.expand(request())
    assert response["cache_hit"] is False
    assert response["result"]["value"] == "0"


@pytest.mark.parametrize(
    "mode,exception",
    [
        ("nonzero", BackendError),
        ("missing", BackendProtocolError),
        ("malformed", BackendProtocolError),
        ("incomplete", BackendProtocolError),
        ("wrong_order", BackendProtocolError),
        ("error", BackendError),
        ("changed_version", BackendProtocolError),
        ("duplicate", BackendProtocolError),
        ("nonfinite", BackendProtocolError),
        ("boolean_protocol", BackendProtocolError),
    ],
)
def test_transport_failures_are_closed(stub_backend, mode, exception):
    with pytest.raises(exception):
        stub_backend.expand(request(stub_mode=mode))


def test_timeout_is_closed(stub_backend):
    # Direct transport request avoids including the health probe in the short limit.
    with pytest.raises(BackendTimeoutError):
        stub_backend._run({"operation": "info", "stub_mode": "timeout"})


@pytest.mark.parametrize("value", [True, 0, -1, float("nan"), float("inf"), "10"])
def test_invalid_timeout(value):
    with pytest.raises(ValueError):
        BackendConfig(timeout=value)


def test_environment_configuration(monkeypatch, tmp_path):
    monkeypatch.setenv("HYPRA_WOLFRAM_EXECUTABLE", "/path/with spaces/kernel")
    monkeypatch.setenv("HYPRA_HYPEXP_PATH", "/packages/HypExp.m")
    monkeypatch.setenv("HYPRA_TIMEOUT", "30")
    monkeypatch.setenv("HYPRA_CACHE_DIR", str(tmp_path))
    monkeypatch.setenv("HYPRA_CACHE", "0")
    config = BackendConfig.from_env()
    assert config.executable == "/path/with spaces/kernel"
    assert config.timeout == 30
    assert config.resolved_cache_dir() == tmp_path
    assert config.use_cache is False


def test_bridge_does_not_evaluate_source_strings():
    bridge = Path(__file__).parents[1] / "src/hypra/wolfram_bridge.wls"
    source = bridge.read_text()
    assert "ToExpression[" not in source
    assert 'Unique["HypraInput`v"]' in source
    assert "HypergeometricPFQ" in source


def test_bridge_delimiters_balance():
    """A static typo guard, explicitly NOT Wolfram execution/validation."""
    source = (Path(__file__).parents[1] / "src/hypra/wolfram_bridge.wls").read_text()
    stack = []
    index, comments = 0, 0
    quoted = False
    while index < len(source):
        char, pair = source[index], source[index : index + 2]
        if comments:
            if pair == "(*":
                comments += 1
                index += 2
                continue
            if pair == "*)":
                comments -= 1
                index += 2
                continue
        elif quoted:
            if char == "\\":
                index += 2
                continue
            if char == '"':
                quoted = False
        elif pair == "(*":
            comments = 1
            index += 2
            continue
        elif char == '"':
            quoted = True
        elif char in "([{":
            stack.append(char)
        elif char in ")]}":
            assert stack and "([{".index(stack.pop()) == ")]}".index(char)
        index += 1
    assert not stack and not comments and not quoted


@pytest.mark.parametrize(
    "field,value",
    [
        ("achieved_through", False),
        ("achieved_through", 0.0),
        ("remainder_order", True),
        ("remainder_order", 1.0),
    ],
)
def test_order_metadata_rejects_boolean_and_float(field, value):
    response = {
        "protocol": 1,
        "status": "ok",
        "complete": True,
        "achieved_through": 0,
        "remainder_order": 1,
        "result": {"head": "Integer", "value": "1"},
        "order_basis": "hypra-pole-bound-plus-hypexp-series",
    }
    response[field] = value
    with pytest.raises(BackendProtocolError):
        WolframBackend._validate_expansion(response, {"through": 0})


def test_cache_rejects_duplicate_keys(stub_backend):
    stub_backend.expand(request())
    cache = next(stub_backend.config.resolved_cache_dir().glob("*.json"))
    text = cache.read_text()
    cache.write_text(text.replace('"key":', '"key":"duplicate", "key":', 1))
    assert stub_backend.expand(request())["cache_hit"] is False


def test_cache_rejects_nonfinite_values(stub_backend):
    stub_backend.expand(request())
    cache = next(stub_backend.config.resolved_cache_dir().glob("*.json"))
    text = cache.read_text()
    cache.write_text(text[:-1] + ', "extra": NaN}')
    assert stub_backend.expand(request())["cache_hit"] is False


def test_kernel_init_fingerprint_covers_sibling_sources(tmp_path):
    kernel = tmp_path / "HypExp" / "Kernel"
    kernel.mkdir(parents=True)
    init = kernel / "init.m"
    init.write_text("(* loader *)")
    source = kernel.parent / "HypExp.m"
    source.write_text("(* version 1 *)")
    before = _package_fingerprint(str(init))
    source.write_text("(* version 2 *)")
    assert _package_fingerprint(str(init)) != before


def test_constant_output_does_not_bypass_hypexp_remainder_guard():
    """Static guard pending execution with a real Wolfram backend."""
    source = (Path(__file__).parents[1] / "src/hypra/wolfram_bridge.wls").read_text()
    assert "requireSeriesDepth[converted, eps, internal, !FreeQ[h, eps]]" in source
    assert "parts === {} && TrueQ[mustHaveRemainder], fail[" in source
