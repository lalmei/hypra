"""Configuration for an externally installed, licensed Wolfram/HypExp backend."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class BackendConfig:
    """Settings for the external Wolfram comparison backend.

    Package paths identify trusted local source files (or package directories).
    These packages execute in the Wolfram kernel: do not configure untrusted files.

    Attributes:
        executable: Path to the Wolfram executable. Omit for auto-discovery.
        hypexp_path: Local HypExp package file or directory.
        hpl_path: Local HPL package file or directory.
        timeout: Kernel subprocess timeout in seconds.
        cache_dir: Override the default external cache directory.
        use_cache: Reuse results from the local disk cache.
    """

    executable: str | None = None
    hypexp_path: str | None = None
    hpl_path: str | None = None
    timeout: float = 120.0
    cache_dir: str | None = None
    use_cache: bool = True

    def __post_init__(self) -> None:
        if (
            isinstance(self.timeout, bool)
            or not isinstance(self.timeout, (int, float))
            or not 0 < self.timeout <= 86400
        ):
            raise ValueError("timeout must be a positive number of seconds, at most 86400")
        for name in ("executable", "hypexp_path", "hpl_path", "cache_dir"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value or "\x00" in value):
                raise ValueError(f"{name} must be a nonempty path string or None")

    @classmethod
    def from_env(cls) -> BackendConfig:
        """Load the external backend settings from Hypra environment variables."""
        return cls(
            executable=os.environ.get("HYPRA_WOLFRAM_EXECUTABLE") or None,
            hypexp_path=os.environ.get("HYPRA_HYPEXP_PATH") or None,
            hpl_path=os.environ.get("HYPRA_HPL_PATH") or None,
            timeout=float(os.environ.get("HYPRA_TIMEOUT", "120")),
            cache_dir=os.environ.get("HYPRA_CACHE_DIR") or None,
            use_cache=os.environ.get("HYPRA_CACHE", "1") not in {"0", "false", "no"},
        )

    def resolved_executable(self) -> str | None:
        """Find the configured or auto-discovered Wolfram executable."""
        if self.executable:
            candidate = os.path.expanduser(self.executable)
            found = shutil.which(candidate)
            return str(Path(found).resolve()) if found else None
        for name in ("wolframscript", "WolframKernel", "MathKernel"):
            if found := shutil.which(name):
                return str(Path(found).resolve())
        candidates = (
            "/Applications/Mathematica.app/Contents/MacOS/WolframKernel",
            "/Applications/Wolfram Mathematica.app/Contents/MacOS/WolframKernel",
            "/Applications/Wolfram.app/Contents/MacOS/WolframKernel",
            "/usr/local/Wolfram/Mathematica/Executables/WolframKernel",
        )
        for candidate in candidates:
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return candidate
        return None

    def resolved_cache_dir(self) -> Path:
        """Return the cache location after expanding user paths."""
        if self.cache_dir:
            return Path(self.cache_dir).expanduser()
        return Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache"))) / "hypra"
