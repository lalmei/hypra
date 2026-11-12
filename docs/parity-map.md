# Native parity map

This maps the public mathematical operations advertised by the upstream HPL
package to Python-native interfaces. It is a capability map, not a claim of
byte-for-byte output compatibility or exhaustive upstream regression equivalence.
The reference source is read only; no upstream reduction tables are copied.

| Mathematical operation | Python implementation | Current boundary |
|---|---|---|
| Hypergeometric regulator expansion | `expand_epsilon` | HypExp 2 listed class patterns; exact affine parameters; explicit domain and resource guards |
| HPL objects and derivatives | `HPL`, SymPy differentiation | Expanded letters 0,±1; no eager branch rewrite |
| Compressed index conversion | `expand_indices`, `compress_indices`, `HPL.from_compressed` | Exact signed integers; bounded output weight |
| Plus/minus kernel notation | `hpl_pm`, `pm_to_letters`, `letters_to_pm` | Exact multilinear conversion; bounded output size |
| Endpoint constant products | `stuffle_expand_constants` | Exact strict colored stuffle algebra, total weight ≤12; not divergent regularization |
| Product expansion | `shuffle_expand`, `shuffle_product` | Same-argument shuffle algebra |
| Logarithm extraction | `extract_hpl_logs` | Leading 1 and trailing 0; retains branch data |
| Shuffle basis reduction | `lyndon_reduce` | Table-free Lyndon basis through weight eight, not upstream's particular chosen basis |
| Nielsen notation | `NielsenPolylog`, `hpl_to_nielsen`, `nielsen_to_hpl` | Exact 0^n1^p words; no general Nielsen-basis claim |
| Elementary/classical conversion | `reduce_known_hpl`; `classical_hpl` | Standard all-weight simple words; all weight-two and weight-three words with at most two distinct letters; other mixed words retained |
| Argument transformations | Negation, reflection, Möbius, inversion, square, fractional maps | Explicit alphabet/domain/boundary/weight contracts |
| Expression argument unification | `transform_hpl_arguments` | Dispatches supported maps; unknown arguments rejected |
| Rational-HPL integration | `primitive_hpl`, `integrate_from_zero` | Poles only 0,±1; definite origin limit must exist |
| Series around zero | `hpl_series`, supported SymPy `.series()` | Explicit logarithmic local-series conventions |
| Values at one | `hpl_at_one`, `hpl_alternating_at_one` | Convergent values only; exact MZV/EulerSum remainder constants |
| Explicit shuffle regularization | `hpl_endpoint_polynomial`, `hpl_expression_endpoint_polynomial`, `shuffle_regularized_hpl_at_one` | Real endpoint polynomial in T=−log(1−z); not an ordinary finite limit |
| Numerical endpoint constants | `evaluate_mzv`, `evaluate_euler_sum`, `.evalf()` | Through weight six; guarded precision evidence, not interval certification |
| Values at i | `hpl_at_i`; numerical `continue_hpl` | Exact all ternary words through weight two, repeated letters and simple classical words through weight 24; higher mixed words are not covered by a special-value table |
| Numerical HPL evaluation | `hpl_evalf`, `evaluate_hpl`, `continue_hpl` | Real interval, explicit slit disk, or explicit off-cut radial path |
| Numerical real-cut boundary values | `evaluate_hpl_boundary` | Exact real argument plus explicit region and upper/lower side; existing map/weight bounds; no guessed small imaginary offset |

## Remaining distinctions

Generic mixed words are not all converted into classical
functions. Named higher-weight constant tables and unrestricted winding paths
are not imported. Numerical stability checks do not certify floating-point
roundoff. Resource limits are deliberate failures, not silent truncations.
The exact hypergeometric domain and numerical HPL domain are separate contracts.
See the capability map and each API's derivation for the precise restrictions.

## HPL 2.0-specific audit

The pinned standalone HPL 2.0 archive is dated **03-November-2011**, not the
older HPL copy bundled with HypExp. The complex-extension paper,
[Maître, hep-ph/0703052](https://arxiv.org/abs/hep-ph/0703052), explicitly
extends HPL to complex arguments. Its feature surface must not be reduced to
"supports harmonic polylogarithms" or confused with hypergeometric expansion.
The [source audit](upstream-source.md#hpl-20-functionality-audit) records the
public declarations and optional table loaders inspected without execution.

The following distinctions remain important even after the native at-i and
real-cut additions:

- **Complex regions and paths:** native argument maps and real-region
  continuation require explicit domains and cut sides. Radial numerical
  continuation follows one specified path. This is not unrestricted complex
  winding/monodromy or automatic branch simplification.
- **Values at i:** an unevaluated `HPL(word, I)` is a lossless symbolic object,
  not an evaluated constant. `hpl_at_i` provides proved low-weight and simple
  all-weight reductions, sometimes retaining exact classical dilogarithms.
  It does not reproduce upstream's weight-dependent special-value tables.
- **Constant reduction:** exact `MZV`/`EulerSum` remainders, shuffle/stuffle
  identities and guarded numerical evaluation are useful native operations,
  but are not a complete reduction to the upstream named constant basis.
- **Minimal sets and known functions:** the native Lyndon basis is a chosen
  table-free shuffle basis, not upstream's particular minimal-set tables.
  Classical conversion still leaves genuinely ternary weight-three words
  and general higher-weight words intact. An exact Nielsen-polylogarithm
  notation/conversion layer covers words 0^n 1^p; it does not reduce arbitrary
  HPL words to Nielsen functions.
- **Numerical assurance:** cut-boundary evaluation maps the exact argument
  before numerics. Its component tail bounds are not a propagated final-error
  certificate; higher-precision stability is checked separately. Neither this
  nor radial continuation provides interval certification of all roundoff.
- **API behavior:** Hypra deliberately avoids HPL 2's global automatic-rewrite
  switches and default cut-side choices. Explicit Python operations are not
  Mathematica syntax/API compatibility. No real upstream-kernel comparison
  suite has been run.

### Next substantive algorithm

The smallest useful remaining symbolic milestone is **all mixed-sign ternary
weight-three HPLs in classical functions on the real interval**. Over the
ordered alphabet `-1 < 0 < 1`, weight three has eight Lyndon basis words. Two
genuinely ternary representatives remain beyond the implemented two-letter primitives:
`(-1,0,1)` and `(-1,1,0)`.

A native implementation can derive the remaining two logarithm/dilogarithm/trilogarithm
primitives, fix each constant by the origin limit, and reconstruct all 27
words through the existing shuffle algebra. Acceptance should check every
exact differentiated recursion, the origin normalization and independent
numerical values. Branch-safe real argument choices are part of the proof;
blind substitution across polylogarithm cuts is not acceptable. This is a
proposed next algorithm, **not an implemented feature claim**.
