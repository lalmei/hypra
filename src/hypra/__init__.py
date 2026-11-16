"""Hypra: a SymPy-native hypergeometric expansion project.

Native exact regulator expansions inspired by HypExp; no Mathematica required.
An experimental external comparison bridge is optional.
"""

from collections.abc import Iterable
from typing import cast

from sympy import hyper

from ._typing import ExactScalar
from .api import ExpansionResult, backend_info, expand_epsilon
from .backend import (
    BackendError,
    BackendProtocolError,
    BackendTimeoutError,
    BackendUnavailableError,
)
from .config import BackendConfig
from .errors import HypraError, IncompleteExpansionError, UnsupportedExpressionError
from .expansion_library import ExpansionLibrary
from .hpl.algebra import hpl_series, integrate_from_zero, primitive_hpl, shuffle_expand
from .hpl.alternating import hpl_alternating_at_one, reduce_euler_sum
from .hpl.arguments import transform_hpl_arguments
from .hpl.boundary_numeric import HPLBoundaryNumericResult, evaluate_hpl_boundary
from .hpl.classical import classical_hpl
from .hpl.constants import ConstantNumericResult, evaluate_euler_sum, evaluate_mzv
from .hpl.endpoints import hpl_at_one, reduce_mzv
from .hpl.inversion import invert_hpl
from .hpl.limits import hpl_limit_at_one
from .hpl.maps import fractional_hpl, square_hpl
from .hpl.mobius import mobius_hpl
from .hpl.named_constants import HPLs6, HPLs7a, HPLs7b, HPLs8a, HPLs8b, HPLs8c, HPLs8d
from .hpl.notation import compress_indices, expand_indices, hpl_pm, letters_to_pm, pm_to_letters
from .hpl.numeric import (
    HPLNumericalConvergenceError,
    HPLNumericResult,
    evaluate_hpl,
    hpl_evalf,
)
from .hpl.radial_continuation import HPLContinuationResult, continue_hpl
from .hpl.real_continuation import analytic_continue_hpl
from .hpl.reduction import extract_hpl_logs, lyndon_factors, lyndon_reduce
from .hpl.regularization import (
    hpl_endpoint_polynomial,
    hpl_expression_endpoint_polynomial,
    shuffle_regularized_hpl_at_one,
)
from .hpl.series_one import hpl_series_at_one
from .hpl.special import hpl_at_i
from .hpl.stuffle import stuffle_expand_constants
from .hpl.taylor import hpl_taylor
from .hpl.transforms import negate_hpl, reflect_hpl
from .hypexp_integrals import beta_log_moment, hypergeometric_log_moment
from .native import NativeBackend, NativeConfig
from .nielsen import NielsenPolylog, hpl_to_nielsen, nielsen_to_hpl, reduce_nielsen
from .serialization import SerializationError, decode_expr, encode_expr
from .symbolic import HPL, MZV, EulerSum

__version__ = "0.2.0"


def hypergeometric(
    upper: Iterable[ExactScalar],
    lower: Iterable[ExactScalar],
    argument: ExactScalar,
) -> hyper:
    """Construct an unevaluated generalized hypergeometric function.

    Construction preserves matching upper and lower parameters for later
    validation by `expand_epsilon`.

    Args:
        upper: Exact numerator parameters in their supplied order.
        lower: Exact denominator parameters in their supplied order.
        argument: Exact argument of the function.

    Returns:
        A SymPy hypergeometric atom with parameter cancellation disabled.
    """
    return cast(hyper, hyper(tuple(upper), tuple(lower), argument, evaluate=False))


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
    "ExpansionLibrary",
    "ExpansionResult",
    "HPLBoundaryNumericResult",
    "HPLContinuationResult",
    "HPLNumericResult",
    "HPLNumericalConvergenceError",
    "HPLs6",
    "HPLs7a",
    "HPLs7b",
    "HPLs8a",
    "HPLs8b",
    "HPLs8c",
    "HPLs8d",
    "HypraError",
    "IncompleteExpansionError",
    "NativeBackend",
    "NativeConfig",
    "NielsenPolylog",
    "SerializationError",
    "UnsupportedExpressionError",
    "__version__",
    "analytic_continue_hpl",
    "backend_info",
    "beta_log_moment",
    "classical_hpl",
    "compress_indices",
    "continue_hpl",
    "decode_expr",
    "encode_expr",
    "evaluate_euler_sum",
    "evaluate_hpl",
    "evaluate_hpl_boundary",
    "evaluate_mzv",
    "expand_epsilon",
    "expand_indices",
    "extract_hpl_logs",
    "fractional_hpl",
    "hpl_alternating_at_one",
    "hpl_at_i",
    "hpl_at_one",
    "hpl_endpoint_polynomial",
    "hpl_evalf",
    "hpl_expression_endpoint_polynomial",
    "hpl_limit_at_one",
    "hpl_pm",
    "hpl_series",
    "hpl_series_at_one",
    "hpl_taylor",
    "hpl_to_nielsen",
    "hypergeometric",
    "hypergeometric_log_moment",
    "integrate_from_zero",
    "invert_hpl",
    "letters_to_pm",
    "lyndon_factors",
    "lyndon_reduce",
    "mobius_hpl",
    "negate_hpl",
    "nielsen_to_hpl",
    "pm_to_letters",
    "primitive_hpl",
    "reduce_euler_sum",
    "reduce_mzv",
    "reduce_nielsen",
    "reflect_hpl",
    "shuffle_expand",
    "shuffle_regularized_hpl_at_one",
    "square_hpl",
    "stuffle_expand_constants",
    "transform_hpl_arguments",
]
