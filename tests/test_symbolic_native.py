"""Exact identities for the native symbolic HPL representation."""

import pytest
import sympy as s
from sympy.core.function import ArgumentIndexError

from hypra.symbolic import HPL, reduce_known_hpl

Z, X = s.symbols("z x")


@pytest.mark.parametrize("letter,kernel", [(-1, 1 / (1 + Z)), (0, 1 / Z), (1, 1 / (1 - Z))])
def test_first_letter_derivative(letter, kernel):
    assert s.diff(HPL((letter, 0, 1), Z), Z) == kernel * HPL((0, 1), Z)
    assert s.diff(HPL((letter,), Z), Z) == kernel


def test_empty_word_derivative_is_zero():
    assert s.diff(HPL((), Z), Z) == 0


def test_derivative_chain_rule():
    assert s.diff(HPL((0, 1), Z**2), Z) == 2 * HPL((1,), Z**2) / Z
    assert s.diff(HPL((-1,), 1 - Z), Z) == -1 / (2 - Z)
    assert s.diff(HPL((1,), Z), X) == 0


def test_discrete_word_is_not_a_differentiable_argument():
    with pytest.raises(ArgumentIndexError):
        HPL((0, 1), Z).fdiff(1)


@pytest.mark.parametrize("weight", range(1, 6))
def test_repeated_letter_reductions(weight):
    assert reduce_known_hpl(HPL((0,) * weight, Z)) == s.log(Z)**weight / s.factorial(weight)
    assert reduce_known_hpl(HPL((1,) * weight, Z)) == (-s.log(1-Z))**weight / s.factorial(weight)
    assert reduce_known_hpl(HPL((-1,) * weight, Z)) == s.log(1+Z)**weight / s.factorial(weight)


@pytest.mark.parametrize("weight", range(2, 6))
def test_classical_polylogarithm_reduction(weight):
    assert reduce_known_hpl(HPL((0,) * (weight - 1) + (1,), Z)) == s.polylog(weight, Z)


def test_empty_and_unknown_word_reduction():
    unknown = HPL((1, 0, -1), Z)
    assert reduce_known_hpl(HPL((), Z)) == 1
    assert reduce_known_hpl(unknown + HPL((0, 1), Z)) == unknown + s.polylog(2, Z)


@pytest.mark.parametrize("word", [(0,), (1,), (-1,), (0, 0), (1, 1), (-1, -1), (0, 1)])
def test_reduction_agrees_with_differential_equation(word):
    direct = s.diff(reduce_known_hpl(HPL(word, Z)), Z)
    recursive = reduce_known_hpl(s.diff(HPL(word, Z), Z))
    assert s.simplify(s.expand_func(direct - recursive)) == 0


def test_reduction_is_explicit_and_preserves_argument():
    atom = HPL((0, 1), 1 - Z)
    assert isinstance(atom, HPL)
    assert reduce_known_hpl(atom) == s.polylog(2, 1-Z)
    assert reduce_known_hpl(3) == 3
    with pytest.raises(TypeError):
        reduce_known_hpl("HPL((0,1),z)")
