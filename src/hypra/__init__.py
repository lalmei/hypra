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
from .hpl_algebra import hpl_series, integrate_from_zero, primitive_hpl, shuffle_expand
from .hpl_alternating import hpl_alternating_at_one, reduce_euler_sum
from .hpl_arguments import transform_hpl_arguments
from .hpl_constants import ConstantNumericResult, evaluate_euler_sum, evaluate_mzv
from .hpl_continue import HPLContinuationResult, continue_hpl
from .hpl_endpoints import hpl_at_one, reduce_mzv
from .hpl_inversion import invert_hpl
from .hpl_maps import fractional_hpl, square_hpl
from .hpl_mobius import mobius_hpl
from .hpl_notation import compress_indices, expand_indices, hpl_pm, letters_to_pm, pm_to_letters
from .hpl_numeric import (
    HPLNumericalConvergenceError,
    HPLNumericResult,
    evaluate_hpl,
    hpl_evalf,
)
from .hpl_reduce import extract_hpl_logs, lyndon_factors, lyndon_reduce
from .hpl_transform import negate_hpl, reflect_hpl
from .native import NativeBackend, NativeConfig
from .serialization import SerializationError, decode_expr, encode_expr
from .symbolic import HPL, MZV, EulerSum

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
    "ConstantNumericResult",
    "EulerSum",
    "ExpansionResult",
    "HPLContinuationResult",
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
    "compress_indices",
    "continue_hpl",
    "decode_expr",
    "encode_expr",
    "evaluate_euler_sum",
    "evaluate_hpl",
    "evaluate_mzv",
    "expand_epsilon",
    "expand_indices",
    "extract_hpl_logs",
    "fractional_hpl",
    "hpl_alternating_at_one",
    "hpl_at_one",
    "hpl_evalf",
    "hpl_pm",
    "hpl_series",
    "hypergeometric",
    "integrate_from_zero",
    "invert_hpl",
    "letters_to_pm",
    "lyndon_factors",
    "lyndon_reduce",
    "mobius_hpl",
    "negate_hpl",
    "pm_to_letters",
    "primitive_hpl",
    "reduce_euler_sum",
    "reduce_mzv",
    "reflect_hpl",
    "shuffle_expand",
    "square_hpl",
    "transform_hpl_arguments",
]
