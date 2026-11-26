import pytest
import sympy as s

from hypra import HPL, HPLPM, hpl_at_minus_one

pytestmark = [pytest.mark.real_backend, pytest.mark.usefixtures("require_real_backend")]


def test_scalars_and_negative_endpoints_against_hpl_2(hpl_reference):
    commutators, endpoints = hpl_reference(
        r"""
HPL`$HPLAutoConvertToKnownFunctions = False;
HPL`$HPLAutoReduceToMinimalSet = False;
HPL`$HPLAutoLogExtract = False;
HPL`$HPLAutoProductExpand = False;
a = HPL`HPL[{1,-1}, x]; b = HPL`HPL[{2}, x];
p = HPL`HPL[{HPL`plus, HPL`minus}, x];
testWords = {{0,1,-1}, {0,-1}, {1,0}, {0,1,0}, {1,0,0}};
pack[v_] := Module[{n = N[v,55]},
  If[!NumberQ[n], Exit[2]];
  ToString[Round[# * 10^42], InputForm] & /@ {Re[n], Im[n]}
];
endpoint[word_, side_] := Module[{v = HPL`HPL[HPL`HPLAtoM[word], x]},
  HPL`HPLAnalyticContinuation[v, HPL`AnalyticContinuationRegion -> HPL`m1to0,
      HPL`AnalyticContinuationSign -> side] /. x -> -1
];
Export[Last[argv], {{a*b-b*a, a*p-p*a},
  Table[{word, pack[endpoint[word,1]], pack[endpoint[word,-1]]}, {word,testWords}]},
  "RawJSON"];
Exit[0];
""",
    )
    x = s.Symbol("x")
    a, b, p = HPL((1, -1), x), HPL((0, 1), x), HPLPM(("plus", "minus"), x)
    assert commutators == [a * b - b * a, a * p - p * a] == [0, 0]
    for word, upper, lower in endpoints:
        for side, (real, imag) in (("upper", upper), ("lower", lower)):
            expected = (s.Integer(real) + s.I * s.Integer(imag)) / 10**42
            actual = hpl_at_minus_one(word, boundary=side).evalf(50)
            assert abs(actual - expected) < s.Rational(1, 10**42)
