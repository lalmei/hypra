"""Hypra: a SymPy-native hypergeometric expansion project.

Native exact regulator expansions inspired by HypExp; no Mathematica required.
An experimental external comparison bridge is optional.
"""

from sympy import hyper

from .api import ExpansionResult, backend_info, expand_epsilon
from .backend import (
    BackendError,
    BackendProtocolError,
    BackendTimeoutError,
    BackendUnavailableError,
)
from .config import BackendConfig
from .errors import HypraError, IncompleteExpansionError, UnsupportedExpressionError
from .hpl_endpoints import hpl_at_one, reduce_mzv
from .native import NativeBackend, NativeConfig
from .serialization import SerializationError, decode_expr, encode_expr
from .symbolic import HPL, MZV

__version__ = "0.1.0a0"


def hypergeometric(upper, lower, argument):
    """Construct an unevaluated generalized hypergeometric function.

    Parameters and the argument use SymPy's exact symbolic representation.
    This constructor does not perform an epsilon expansion.
    """
    return hyper(tuple(upper), tuple(lower), argument, evaluate=False)


__all__ = [
    "HPL",
    "MZV",
    "BackendConfig",
    "BackendError",
    "BackendProtocolError",
    "BackendTimeoutError",
    "BackendUnavailableError",
    "ExpansionResult",
    "HypraError",
    "IncompleteExpansionError",
    "NativeBackend",
    "NativeConfig",
    "SerializationError",
    "UnsupportedExpressionError",
    "__version__",
    "backend_info",
    "decode_expr",
    "encode_expr",
    "expand_epsilon",
    "hpl_at_one",
    "hypergeometric",
    "reduce_mzv",
]
