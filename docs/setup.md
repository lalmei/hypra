# Native Python setup

```sh
cd ~/projects/hypra
uv sync
uv run hypra backend-info
uv run pytest -q
uv run python scripts/validate.py --output-dir reports/local
```

No Mathematica, Wolfram license or API key is required. Native Python/SymPy is
always the default and unsupported native cases never fall back to another
backend. `backend_info()` reports the native engine and its version; each
expansion records exact input, assumptions, requested/internal order, source
fingerprint, SymPy version and cache status.

The Python quality workflow runs alongside the test matrix. Run the same
checks locally from the repository root.

```sh
uv run --locked --group quality ruff format --check src tests examples scripts
uv run --locked --group quality ruff check .
uv run --locked --group quality pyright --pythonpath .venv/bin/python
```

```python
from hypra import NativeConfig, expand_epsilon
config = NativeConfig(max_internal_order=12, use_cache=True,
                      cache_dir="/your/cache/path")
# expand_epsilon(expr, eps, through=3, assumptions=domain, config=config)
```

Default native cache: `~/.cache/hypra/native`. Keys include the complete exact
request, assumptions, orders and algorithm/source/SymPy versions. Checksums
reject accidental corruption; this is not authentication against a malicious
local user controlling the cache. `use_cache=False` recomputes everything.

## Optional external comparison bridge

The optional Wolfram adapter runs trusted local HypExp and HPL packages
under a licensed Wolfram kernel. Select `backend="wolfram"` explicitly.
The live tests exercise selected cases on Wolfram Engine 15.

- `HYPRA_WOLFRAM_EXECUTABLE` sets the executable path. For Wolfram Engine
  on macOS, use `/Applications/Wolfram Engine.app/Contents/MacOS/WolframKernel`.
- `HYPRA_HYPEXP_PATH` and `HYPRA_HPL_PATH` identify trusted local package
  entry points.
- `HYPRA_TIMEOUT` sets the process deadline in seconds, default 120.
- `HYPRA_CACHE_DIR` and `HYPRA_CACHE=0` control the external result cache.

`uv run hypra backend-info --backend wolfram` checks the configured installation.
To run the live suite with an activated Wolfram Engine on macOS:

```sh
HYPRA_WOLFRAM_EXECUTABLE='/Applications/Wolfram Engine.app/Contents/MacOS/WolframKernel' \
HYPRA_RUN_REAL=1 uv run pytest -q -m real_backend
```

The suite disables caching. If both package paths are unset, it uses the
extracted archives under `reference/upstream` when available. It creates
HypExp's initial library index in a temporary archive copy. Missing kernel
or package setup raises a test error.

## Native modules and verification

No additional system software is needed for integer reductions, matched-half
colored sums, canonical half-integer ODE systems, rational-HPL integration or
endpoint rewriting. They are Python modules in the same package. Their exact
contracts and remaining gaps are listed in [scope](scope.md); roadmap items
are not silently treated as supported requests.

Unit tests live in `tests/unit/`, grouped by package for HPL, native, and
symbolic code. Integration tests live in `tests/integration/` and exercise the
CLI, public expansions, validation examples, and the optional Wolfram bridge.
Run a group directly with pytest.

```sh
uv run pytest -q tests/unit
uv run pytest -q tests/unit/hpl
uv run pytest -q tests/integration -m 'not real_backend'
```

For focused offline mathematical regressions, run the following tests.

```sh
uv run pytest -q \
  tests/unit/native/test_integer_shift_math.py tests/unit/native/test_half_system.py \
  tests/unit/native/test_endpoints.py tests/unit/hpl/test_hpl_algebra.py \
  tests/unit/native/test_nested_sums.py
```

The broader test suite should still be run before accepting a change: kernel
checks alone do not establish public dispatch, cache invalidation or complete
Laurent-prefactor behavior. Compiling documentation LaTeX to PDF is optional
and separate from running the package; `result.latex()` needs no TeX install.

Native resource guards bound regulator depth, shift count and selected algebraic
sizes. For an execution deadline, use `NativeConfig(timeout_seconds=30)`.
This runs the exact JSON request in a separate Python process and terminates it
on timeout; no partial expansion is returned. Process startup counts toward the
deadline. Without this option, native calculations run in-process without a
wall-clock deadline. Native and optional external-adapter timeouts are separate.
The deadline covers the native backend calculation, not earlier input validation,
assumption checks or endpoint preprocessing in the public API.
