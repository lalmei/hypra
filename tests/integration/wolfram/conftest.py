"""Set up the pinned packages for tests that run the real Wolfram kernel."""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from hypra import BackendConfig, backend_info


@pytest.fixture(scope="module")
def require_real_backend(tmp_path_factory):
    if os.environ.get("HYPRA_RUN_REAL") != "1":
        pytest.skip("Real backend is opt-in: HYPRA_RUN_REAL=1")
    with pytest.MonkeyPatch.context() as env:
        env.setenv("HYPRA_CACHE", "0")
        if not (os.environ.get("HYPRA_HYPEXP_PATH") or os.environ.get("HYPRA_HPL_PATH")):
            upstream = Path(__file__).resolve().parents[3] / "reference" / "upstream"
            hpl = upstream / "HPL-2.0" / "HPL.m"
            hypexp = upstream / "HypExp-2.0"
            if hpl.is_file() and (hypexp / "HypExp.m").is_file():
                workspace = tmp_path_factory.mktemp("upstream-hypexp") / "HypExp-2.0"
                shutil.copytree(hypexp, workspace)
                known = workspace / ("known.m" if os.name == "nt" else ".known.m")
                if not known.exists():
                    known.write_text("{KnownToOrder1}\n", encoding="utf-8")
                env.setenv("HYPRA_HPL_PATH", str(hpl))
                env.setenv("HYPRA_HYPEXP_PATH", str(workspace / "HypExp.m"))
        info = backend_info(backend="wolfram")
        assert info["available"], info["error"]
        yield info


@pytest.fixture
def hpl_reference(tmp_path, require_real_backend):
    """Run a comparison directly against the configured HPL package."""
    config = BackendConfig.from_env()
    executable = config.resolved_executable()
    assert executable is not None
    script = tmp_path / "compare_hpl.wls"
    output = tmp_path / "values.json"

    def compare(source):
        script.write_text(
            r"""argv = If[$ScriptCommandLine === {}, $CommandLine, $ScriptCommandLine];
Global`$HPLPath = DirectoryName[argv[[-2]]];
Get[argv[[-2]]];
words = Flatten[Table[Tuples[{-1, 0, 1}, n], {n, 1, 6}], 1];
toPM[word_] := Replace[word, {-1 -> HPL`minus, 1 -> HPL`plus}, {1}];
fromPM[word_] := Replace[word, {HPL`minus -> -1, HPL`plus -> 1}, {1}];
"""
            + source,
            encoding="utf-8",
        )
        mode = "-file" if "wolframscript" in Path(executable).name.lower() else "-script"
        run = subprocess.run(
            [executable, mode, str(script), require_real_backend["hpl_path"], str(output)],
            capture_output=True,
            text=True,
            timeout=max(config.timeout, 180),
            check=False,
        )
        assert run.returncode == 0 and output.is_file(), run.stdout + run.stderr
        return json.loads(output.read_text())

    return compare
