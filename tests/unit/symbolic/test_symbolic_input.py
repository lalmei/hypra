from sympy import Rational, Symbol, hyper

from hypra import hypergeometric


def test_preserves_epsilon_parameters():
    eps, z = Symbol("eps"), Symbol("z")
    result = hypergeometric([1 + eps, 1], [2 - eps], z)
    assert result == hyper((1 + eps, 1), (2 - eps,), z, evaluate=False)
    assert eps in result.free_symbols


def test_preserves_exact_half_integer_parameter():
    z = Symbol("z")
    result = hypergeometric([Rational(1, 2)], [1], z)
    assert result.ap[0] == Rational(1, 2)


def test_does_not_cancel_parameters_at_construction():
    z = Symbol("z")
    result = hypergeometric([1], [1], z)
    assert result.ap == (1,)
    assert result.bq == (1,)
