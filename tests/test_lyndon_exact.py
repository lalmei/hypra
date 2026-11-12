"""Formal shuffle-algebra checks, independent of HPL argument series."""

from itertools import pairwise, product

import pytest
import sympy as s

from hypra.hpl_algebra import shuffle_expand
from hypra.hpl_reduce import lyndon_factors, lyndon_reduce
from hypra.symbolic import HPL

Z = s.Symbol("z")
WORDS = list(product((-1, 0, 1), repeat=4)) + [
    (1, 0, -1, 0, 1),
    (1, 1, 0, -1, -1),
    (0, -1, 0, -1, 0),
    (-1, 0, -1, 0, 1),
    (1, 0, 1, 0, 1),
    (1, 1, 1, 1, 1),
    (0, 0, 0, 0, 0),
    (-1, -1, -1, -1, -1),
    (1, -1, 1, -1, 1),
]


@pytest.mark.parametrize("word", WORDS)
def test_lyndon_reduction_exact_shuffle_identity(word):
    factors = lyndon_factors(word)
    assert sum(factors, ()) == word
    assert all(left >= right for left, right in pairwise(factors))
    # Independent defining Lyndon criterion: strictly below each proper suffix.
    for factor in factors:
        assert all(factor < factor[k:] for k in range(1, len(factor)))
    reduced = lyndon_reduce(HPL(word, Z))
    assert shuffle_expand(reduced, Z) == HPL(word, Z)
    assert lyndon_reduce(reduced) == reduced
    for atom in reduced.atoms(HPL):
        atom_word = tuple(atom.word)
        assert all(atom_word < atom_word[k:] for k in range(1, len(atom_word)))


def test_equal_lyndon_factor_multiplicity():
    # Two identical nontrivial factors produce a factorial diagonal, not one.
    word = (-1, 0, -1, 0)
    expected = HPL((-1, 0), Z) ** 2 / 2 - 2 * HPL((-1, -1, 0, 0), Z)
    assert s.expand(lyndon_reduce(HPL(word, Z)) - expected) == 0
