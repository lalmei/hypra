"""Exact independent closed-form checks, without Mathematica or numeric fitting."""
import importlib.util
import sys
from pathlib import Path

import pytest
import sympy as s

from hypra import NativeConfig, expand_epsilon

_PATH=Path(__file__).parents[1]/'examples'/'ported_examples.py'
_SPEC=importlib.util.spec_from_file_location('ported_examples',_PATH)
examples=importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name]=examples
_SPEC.loader.exec_module(examples)


@pytest.mark.parametrize('case',examples.examples(),ids=lambda case:case.name)
def test_exact_independently_derived_examples(case):
    result=expand_epsilon(case.expression,examples.EPSILON,through=case.through,
                          assumptions=case.assumptions,config=NativeConfig(use_cache=False))
    # Shuffle extraction is exact; expand_log uses Z's positive assumption,
    # not an unconditional complex logarithm identity.
    reduced=examples.classical_reduction(result.polynomial)
    assert s.simplify(s.expand(reduced-case.expected))==0


def test_generating_formulas_match_exact_defining_coefficients():
    e,z=examples.EPSILON,examples.Z
    for case in (examples.integer_polylog_example(4),examples.matched_half_example()):
        actual=s.series(case.expected,z,0,4).removeO().expand()
        expected=s.S.Zero
        for n in range(4):
            coefficient=s.prod(s.rf(a,n) for a in case.expression.ap)/(s.factorial(n)*s.prod(s.rf(b,n) for b in case.expression.bq))
            expected+=z**n*s.series(coefficient,e,0,case.through+1).removeO()
        assert s.simplify(s.expand(actual-expected))==0
