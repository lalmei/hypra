# Native capability map

Hypra is a Python/SymPy design inspired by the algorithms and mathematics
of HypExp/HypExp 2, not a literal source-code translation. Mathematica is not
required. An implemented mathematical kernel is not automatically a claim of
complete public API coverage or full upstream parity.

## Components and contracts

| Module | Implemented mathematical work | Important restrictions |
| --- | --- | --- |
| `nested_sums.py` | Strict harmonic-sum product algebra, formal exponential, canonical integer hypergeometric coefficients to HPLs | Canonical upper bases 0/1 and lower base 1; exact rational slopes |
| `integer.py` | Integer-base `pF(p-1)` reduction, upper raising/lower shifting, finite-prefix/tail reindexing | Fixed singular lower parameters rejected; moving lower poles require extra orders; high offsets can be expensive |
| `colored_sums.py` | Signed strict sums, colored stuffle products, even-index projection, canonical matched-half coefficients | Equal half-upper/lower counts; alphabet signs ±1; rational slopes; HPL argument `sqrt(z)` |
| `matched_half.py` | Integer-offset reduction around matched half-integer bases, including prefix/tail reindexing | `pF(p-1)` with equal counts of half-integral upper/lower bases; fixed singular lower rejected |
| `contiguous.py` / `canonical.py` | Generic upper/lower contiguous differential operators, Laurent convolution, one-parameter removable-limit regularization | Supported canonical families only; bounded shift count/depth; unresolved singular limits fail closed |
| `half_system.py` | Canonical rationalized epsilon ODE systems A/B/C | Patterns specified below; rank at least two; symbolic slopes allowed internally, not a general public symbolic-slope promise |
| `hpl_algebra.py` | Shuffle linearization, rational-HPL primitives, definite origin integration, local exact series | Rational poles only 0/±1; finite origin constant must be proved; explicit weight/pole/term budgets |
| `endpoints.py` | Finite terminating sums; convergent Gauss `2F1(1)` rewrite | Fixed nonpositive upper termination only; affine Gauss parameters and proven convergence; no divergent endpoint continuation |

The native dispatcher uses integer and matched-half reductions. The canonical half-system kernels cover:

A. All upper bases zero, exactly one lower base 1/2 and the rest 1.
B. Exactly one upper base 1/2, remaining upper bases zero, all lower bases 1.
C. `2F1` upper bases (0,1/2), lower base 1/2.
D. All upper and lower bases 1/2, arbitrary rank at least two, using a rational gauge transformation.

The generic contiguous engine now supplies integer shifts for the supported
canonical patterns, including unmatched-half A/B at arbitrary rank and all-half families. This covers the canonical parameter-count classes listed
by HypExp 2, not every conceivable half-count pattern or full upstream API.
For example mixed `3F2` with two half-upper and one half-lower remains unsupported. See
[half-system derivation](half-system.md) and [integer identities](integer-port.md).
Endpoint rewriting must run before pole planning to include hidden Gamma
poles; see [endpoint conditions](endpoints.md).

## Exactness, branches and failure behavior

The main public nonterminating domain remains real `0<z<1` with
regulator-independent argument and affine exact parameter data. In case A,
the internal HPL variable is imaginary; this is a documented analytic
continuation of the regular origin germ, not permission for arbitrary
complex input branches. Matched-half functions use principal `sqrt(z)`.

The HPL alphabet convention is H1(z)=-log(1-z). Colored negative signs mean
(-1)^n, not (-1)^(n-1). Shuffle and stuffle are different algebras: HPL products
use shuffle; strict finite sums use stuffle including a merged-index term.

Laurent prefactors and intermediate parameter-shift poles require deeper
coefficient orders. Resource bounds, unsupported structures, unproved origin
limits and unproved endpoint convergence raise errors rather than returning
an incomplete expression as a success. There is no external-backend fallback.

## Evidence and remaining scope

Tests include exact finite coefficient comparisons with the defining
Pochhammer series, independent HPL Taylor expansion, operator identities,
round-trip serialization and selected numeric secondary checks. Finite tests
are not proof of every high-order/rank performance claim. The J4 endpoint
exchange proof applies to J4, not arbitrary integrals.

Uncovered parameter patterns, nonterminating ranks other than `pF(p-1)`,
complete special-value tables, general contour continuation and full upstream
HPL simplification are not implemented. The optional Wolfram bridge remains
unverified comparison infrastructure and is not a native dependency.

## Mathematical sources

- Huber and Maitre, HypExp, https://arxiv.org/abs/hep-ph/0507094:
  expansion around integer parameters and nested-sum methods.
- Huber and Maitre, HypExp 2, https://arxiv.org/abs/0708.2443:
  half-integer classes, reduction and rationalized-variable expansions.
- Standard defining series: https://dlmf.nist.gov/16.2.E1 .
- Gauss summation with its convergence condition: https://dlmf.nist.gov/15.4.E20 .

No equation number from the papers is asserted without an independently
verified match. The module derivations document the identities actually used.
No upstream Mathematica source or lookup tables are redistributed.
