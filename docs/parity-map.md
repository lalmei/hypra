# Native parity map

This maps the public mathematical operations advertised by the upstream HPL
package to Python-native interfaces. It is a capability map, not a claim of
byte-for-byte output compatibility or exhaustive upstream regression equivalence.
The reference source is read only; no upstream reduction tables are copied.
The [port audit](port-audit.md) records fixes, independent review, and live checks.

| Mathematical operation | Python implementation | Current boundary |
|---|---|---|
| Hypergeometric regulator expansion | `expand_epsilon` | HypExp 2 listed class patterns, exact affine parameters and [explicit argument domains](scope.md#exactness-branches-and-failure-behavior) |
| Hypergeometric unit endpoints | `expand_epsilon` with `z=1` | [Convergent Gauss and supported higher-rank cases](endpoints.md), plus opt-in meromorphic parameter continuation |
| HPL objects and derivatives | `HPL`, SymPy differentiation | Expanded letters 0,±1; no eager branch rewrite |
| Compressed index conversion | `expand_indices`, `compress_indices`, `HPL.from_compressed` | Exact signed integers; bounded output weight |
| Plus/minus kernel notation | `hpl_pm`, `pm_to_letters`, `letters_to_pm` | Exact multilinear conversion; bounded output size |
| Endpoint constant products | `stuffle_expand_constants` | Exact strict colored stuffle algebra, total weight ≤12; not divergent regularization |
| Product expansion | `shuffle_expand`, `shuffle_product` | Same-argument shuffle algebra |
| Logarithm extraction | `extract_hpl_logs` | Leading 1 and trailing 0; retains branch data |
| HPL 2.0 minimal sets | `reduce_hpl_to_minimal_set`, `HPLPM` | Exact ordinary and plus/minus bases through weight six, generated from shuffle identities |
| Shuffle basis reduction | `lyndon_reduce` | Table-free Lyndon basis through weight eight, not upstream's particular chosen basis |
| Nielsen notation | `NielsenPolylog`, `hpl_to_nielsen`, `nielsen_to_hpl`, `reduce_nielsen` | Exact 0^n1^p words, numerical evaluation, and proved values at 0,±1 |
| Known-function conversion | `reduce_known_hpl`, `classical_hpl`, `reduce_hpl_to_known_functions` | Repeated letters, all ternary cubic words, two-letter quartic words and signed Nielsen families; other mixed words retained as upstream does |
| Argument transformations | Negation, reflection, Möbius, inversion, square, fractional maps | Explicit alphabet/domain/boundary/weight contracts |
| Expression argument unification | `transform_hpl_arguments` | Dispatches supported maps; unknown arguments rejected |
| Rational-HPL integration | `primitive_hpl`, `integrate_from_zero` | Poles only 0,±1; definite origin limit must exist |
| Near-one series and combined limits | `hpl_series_at_one`, `hpl_limit_at_one` | Explicit logarithmic remainders; bounded weight/order/poles; no divergent-limit assignment |
| Series around zero | `hpl_series`, supported SymPy `.series()` | Explicit logarithmic local-series conventions |
| Values at one | `hpl_at_one`, `hpl_alternating_at_one`, `reduce_mzv`, `reduce_euler_sum` | All convergent MZVs through weight eight; colored Euler sums through weight six; proved general families; exact primitive remainders |
| Explicit shuffle regularization | `hpl_endpoint_polynomial`, `hpl_expression_endpoint_polynomial`, `shuffle_regularized_hpl_at_one` | Real endpoint polynomial in T=−log(1−z); not an ordinary finite limit |
| Named HPL constants | `HPLs6`, `HPLs7a`–`HPLs8d` | Exact weak harmonic sums expressed in the strict Euler-sum convention |
| Numerical endpoint constants | `evaluate_mzv`, `evaluate_euler_sum`, `.evalf()` | Through weight eight; guarded precision evidence, not interval certification |
| Values at i | `hpl_at_i`, `hpl_pm_at_i`; numerical `continue_hpl` | All ternary words through weight six in the shipped named-constant/primitive basis; repeated letters and simple classical words through weight 24 |
| Numerical HPL evaluation | `hpl_evalf`, `evaluate_hpl`, `continue_hpl` | Real interval, explicit slit disk or off-cut principal radial path, also via SymPy `.evalf()` |
| Numerical real-cut boundary values | `evaluate_hpl_boundary` | All words through weight 24; exact real target and explicit region/side; symbolic maps or finite-detour Taylor transport |

## HypExp 2.0 public operations

The public declarations in the pinned `HypExp.m` distinguish expansion,
definite logarithmic integrals, and Mathematica library management.

| Upstream operation | Python interface | Current boundary |
|---|---|---|
| `HypExp`, `MyHypExpansion` | `expand_epsilon` with an exact SymPy `hyper` | See [supported families and domains](scope.md#hypexp-2-class-correspondence) and [unit endpoint coverage](endpoints.md). The class list describes kernel coverage. |
| `HypExpU` | `beta_log_moment(n, m, p)` | Convergence requires `p + m >= 0`. Log weight through 32 and `p <= 4096`, using exact Beta/Gamma coefficients. |
| `HypExpInt` | `hypergeometric_log_moment(p, n, m, k, q, z)` | Convergent endpoints, real mixed-log kernels for `z<1`, and explicit principal nonreal arguments through upstream log weight five, subject to power and algebra budgets. |
| `HypExpAddToLib`, `HypExpIsKnownToOrder` | `ExpansionLibrary.add_to_lib`, `ExpansionLibrary.is_known_to_order` | Persistent symbolic templates indexed by sorted integer bases. Opt in with `NativeConfig(expansion_library=...)` for real interior expansions. |
| `$HypExpFailed` | Typed exceptions | Failures do not masquerade as computed coefficients. |
| `$HypExpPolyLogRules` argument relations | `transform_polylog_arguments` | All six Appendix A.1 cross-ratio arguments for classical/Nielsen functions through weight eight, with explicit cut side and an exact common-argument HPL basis. |
| `$HypExpPolyLogRules` named-function reductions | `reduce_polylogs` | Classical and Nielsen functions through weight five at all six cross-ratio arguments reduce to logarithms, classical polylogarithms and canonical `S32(z)`, with explicit cut sides. |

### Registering expansions

```python
from hypra import ExpansionLibrary, NativeConfig

library = ExpansionLibrary("expansions.sqlite")
library.add_to_lib([0, 1], [1], through=2)
assert library.is_known_to_order([1, 0], [1], through=1)
config = NativeConfig(expansion_library=str(library.path))
```

Pass `config` to `expand_epsilon` to specialize the stored template to the
request's slopes and argument. The planner can request a higher order to
account for Laurent poles. If the library lacks that order or the argument
falls outside the real interior or integer-template scope, expansion uses the usual native kernels.
Registration computes the generic expansion with the native order and
algebra budgets. A failed registration preserves existing entries.

## Remaining distinctions

The Python APIs expose explicit operations instead of Mathematica's global
rewrite switches. Generic higher mixed words and arbitrary winding paths
are outside the corresponding shipped
upstream operations. Numerical stability checks do not certify floating-point
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
- **Values at i:** `HPL(word, I).evalf()` evaluates mixed words numerically
  by principal radial continuation. `hpl_at_i` and `hpl_pm_at_i` reduce all
  words through weight six to named constants and the shipped primitive basis.
  Retained `HPLPM(word,I)` constants have exact AST transport and independent
  numerical evaluation.
- **Constant reduction:** exact shuffle/stuffle systems cover all convergent
  MZVs through weight eight and colored sums through weight six. The latter
  retains the primitive `EulerSum((5,1),(-1,-1)) = HPLs6-zeta(6)`.
  Higher general values remain exact unless a proved family applies.
- **Minimal sets and known functions:** the minimal-set reducer reproduces the
  shipped ordinary and plus/minus bases through weight six. The separate
  Lyndon basis remains available through weight eight.
  Known-function conversion covers all 27 ternary weight-three words,
  two-letter quartic words and signed Nielsen families. Higher generic
  words remain intact. Plus/minus expressions are normalized for algebra,
  transformations, integration and series.
- **Numerical assurance:** cut-boundary evaluation maps the exact argument
  before numerics. Its component tail bounds are not a propagated final-error
  certificate; higher-precision stability is checked separately. Neither this
  nor radial continuation provides interval certification of all roundoff.
- **API behavior:** Hypra deliberately avoids HPL 2's global automatic-rewrite
  switches and default cut-side choices. Explicit Python operations are not
  Mathematica syntax/API compatibility. The opt-in live suite compares
  selected HypExp coefficients and HPL values with the upstream packages.

### Further symbolic coverage

All ternary cubic words now reduce through proved square-log, Möbius and
shuffle identities, with real classical arguments on the stated interval.
Near-one power/log expansions and combined ternary rational-pole limits are
also implemented with explicit remainder and resource contracts. Exact at-i
comparison covers all 1,053 plus/minus words of weights four through six;
numerical checks independently cover all 1,093 ordinary words through six.
