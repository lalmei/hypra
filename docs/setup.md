# Native Python setup

```sh
cd ~/project/hypra
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

The original wrapper direction left an experimental Wolfram adapter in the
repository, but it is **not a dependency or a validated feature**. No Wolfram
installation was found, and source inspection showed that common upstream
HypExp paths remove their `O` marker, which this fail-closed adapter rejects.
Do not mistake mock protocol tests for actual Wolfram validation.

Only if independently developing that adapter, select `backend="wolfram"`
explicitly and configure trusted local software:

- `HYPRA_WOLFRAM_EXECUTABLE`: executable path, not a shell command.
- `HYPRA_HYPEXP_PATH`, `HYPRA_HPL_PATH`: trusted local package entry points.
- `HYPRA_TIMEOUT`: process deadline in seconds, default 120.
- `HYPRA_CACHE_DIR`, `HYPRA_CACHE=0`: optional external cache controls.

`uv run hypra backend-info --backend wolfram` checks that optional installation.
`HYPRA_RUN_REAL=1 uv run pytest -m real_backend` requires real execution and
fails rather than skipping missing setup. This path is for future comparison
work, not a prerequisite to use the native package.

## Native modules and verification

No additional system software is needed for integer reductions, matched-half
colored sums, canonical half-integer ODE systems, rational-HPL integration or
endpoint rewriting. They are Python modules in the same package. Their exact
contracts and remaining gaps are listed in [scope](scope.md); roadmap items
are not silently treated as supported requests.

For focused offline mathematical regressions:

```sh
uv run pytest tests/test_integer_shift_math.py tests/test_half_system.py \
  tests/test_endpoints.py tests/test_hpl_algebra.py tests/test_nested_sums.py -q
```

The broader test suite should still be run before accepting a change: kernel
checks alone do not establish public dispatch, cache invalidation or complete
Laurent-prefactor behavior. Compiling documentation LaTeX to PDF is optional
and separate from running the package; `result.latex()` needs no TeX install.

Native resource guards bound regulator depth, shift count and selected algebraic
sizes. They are not a hard wall-clock timeout: a difficult SymPy simplification
can still take a long time. Run long academic jobs in a separately managed
process if an operating-system execution deadline is required. The optional
external adapter's timeout does not apply to native computations.
