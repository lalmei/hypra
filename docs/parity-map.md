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
| Product expansion | `shuffle_expand`, `shuffle_product` | Same-argument shuffle algebra |
| Logarithm extraction | `extract_hpl_logs` | Leading 1 and trailing 0; retains branch data |
| Shuffle basis reduction | `lyndon_reduce` | Table-free Lyndon basis through weight eight, not upstream's particular chosen basis |
| Elementary/classical conversion | `reduce_known_hpl`; `classical_hpl` | Standard all-weight simple words; all weight-two words; higher mixed words retained |
| Argument transformations | Negation, reflection, Möbius, inversion, square, fractional maps | Explicit alphabet/domain/boundary/weight contracts |
| Expression argument unification | `transform_hpl_arguments` | Dispatches supported maps; unknown arguments rejected |
| Rational-HPL integration | `primitive_hpl`, `integrate_from_zero` | Poles only 0,±1; definite origin limit must exist |
| Series around zero | `hpl_series`, supported SymPy `.series()` | Explicit logarithmic local-series conventions |
| Values at one | `hpl_at_one`, `hpl_alternating_at_one` | Convergent values only; exact MZV/EulerSum remainder constants |
| Numerical endpoint constants | `evaluate_mzv`, `evaluate_euler_sum`, `.evalf()` | Through weight six; guarded precision evidence, not interval certification |
| Values at i | Exact `HPL(word,I)`, numerical `continue_hpl` | No copied special-value tables; explicit radial path |
| Numerical HPL evaluation | `hpl_evalf`, `evaluate_hpl`, `continue_hpl` | Real interval, explicit slit disk, or explicit off-cut radial path |

## Remaining distinctions

Generic mixed words above weight two are not all converted into classical
functions. Named higher-weight constant tables and unrestricted winding paths
are not imported. Numerical stability checks do not certify floating-point
roundoff. Resource limits are deliberate failures, not silent truncations.
The exact hypergeometric domain and numerical HPL domain are separate contracts.
See the capability map and each API's derivation for the precise restrictions.
