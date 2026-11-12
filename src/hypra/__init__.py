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
from .hpl_numeric import (
    HPLNumericalConvergenceError,
    HPLNumericResult,
    evaluate_hpl,
    hpl_evalf,
)
from .hpl_reduce import extract_hpl_logs
from .hpl_transform import negate_hpl, reflect_hpl
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
    "HPLNumericResult",
    "HPLNumericalConvergenceError",
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
    "evaluate_hpl",
    "expand_epsilon",
    "extract_hpl_logs",
    "hpl_at_one",
    "hpl_evalf",
    "hypergeometric",
    "negate_hpl",
    "reduce_mzv",
    "reflect_hpl",
]
