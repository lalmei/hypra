"""Public errors: symbolic failures are never returned as successful expansions."""


class HypraError(Exception):
    """Base error for Hypra."""


class UnsupportedExpressionError(HypraError, ValueError):
    """Input cannot be expanded safely by the current implementation."""


class IncompleteExpansionError(HypraError):
    """Backend did not supply a complete, exact Laurent expansion."""
