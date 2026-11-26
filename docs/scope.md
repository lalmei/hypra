# Native capability map

Hypra is a Python/SymPy design inspired by the algorithms and mathematics
of HypExp/HypExp 2, not a literal source-code translation. Mathematica is not
required. An implemented mathematical kernel is not automatically a claim of
complete public API coverage or full upstream parity.

## Components and contracts

| Module | Implemented mathematical work | Important restrictions |
| --- | --- | --- |
| `native/nested_sums.py` | Strict harmonic-sum product algebra, formal exponential, canonical integer hypergeometric coefficients to HPLs | Canonical upper bases 0/1 and lower base 1; exact rational or guarded polynomial symbolic slopes |
| `native/integer.py` | Integer-base `pF(p-1)` reduction, upper raising/lower shifting, finite-prefix/tail reindexing | Fixed singular lower parameters rejected; moving lower poles require extra orders; high offsets can be expensive |
| `native/colored_sums.py` | Signed strict sums, colored stuffle products, even-index projection, canonical matched-half coefficients | Equal half-upper/lower counts; alphabet signs ±1; guarded exact slopes; HPL argument `sqrt(z)` |
| `native/matched_half.py` | Integer-offset reduction around matched half-integer bases, including prefix/tail reindexing | `pF(p-1)` with equal counts of half-integral upper/lower bases; fixed singular lower rejected |
| `native/contiguous.py` / `native/canonical.py` | Generic upper/lower contiguous differential operators, Laurent convolution, one-parameter removable-limit regularization | Supported canonical families only; bounded shift count/depth; unresolved singular limits fail closed |
| `native/half_system.py` | Canonical rationalized epsilon ODE systems A/B/C/D | Patterns specified below; rank at least two; guarded polynomial symbolic slopes; symbolic singular lower slopes need explicit nonzero assumptions |
| `hpl/algebra.py` | Shuffle linearization, rational-HPL primitives, definite origin integration, local exact series | Rational poles only 0/±1; finite origin constant must be proved; explicit weight/pole/term budgets |
| `native/endpoints.py` / `native/unit_limit.py` | Finite terminating sums, convergent unit endpoints, opt-in meromorphic continuation at `z=1` | Termination requires a fixed nonpositive integer upper parameter. Higher-rank endpoints require supported affine integer or half-integer families and bounded HPL limits |
| `hypexp_integrals.py` | Exact definite beta-log and mixed-log moments | Real kernels and explicit convergence, weight and power bounds. See the [integral APIs](python-api.md#definite-logarithmic-integrals) |

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
Endpoint rewriting runs before Laurent planning so Gamma quotients and
excess-raising recurrences can expose additional regulator poles. The default
requires defining-series convergence. Select
`NativeConfig(unit_endpoint="meromorphic")` to continue in the parameters at
fixed `z=1`. See the [endpoint conditions and budgets](endpoints.md).

## HPL 2-inspired native layer

The [HPL parity map](parity-map.md) records numerical continuation, exact
special values, classical reduction, shuffle/stuffle algebra and explicit
regularization coverage. The [Taylor API](hpl-taylor.md) also integrates
ordinary interior-point expansions with SymPy's `.series()` method. These
are Python-native operations, not a redistributed Mathematica package or a
claim that every upstream special-value table has been reproduced.

## Exactness, branches and failure behavior

The default nonterminating domain accepts provably negative real arguments
and `0<z<1`, with regulator-independent argument and affine exact parameter
data. For a symbolic real argument, give assumptions establishing the chosen
interval. In case A on the positive interval,
the internal HPL variable is imaginary; this is a documented analytic
continuation of the regular origin germ, not permission for arbitrary
complex input branches. Matched-half functions use principal `sqrt(z)`.

The HPL alphabet convention is H1(z)=-log(1-z). Colored negative signs mean
(-1)^n, not (-1)^(n-1). Shuffle and stuffle are different algebras: HPL products
use shuffle; strict finite sums use stuffle including a merged-index term.

Laurent prefactors and intermediate parameter-shift poles require deeper
coefficient orders. Unsupported structures, unproved limits or convergence,
and requests exceeding resource bounds raise errors. The native backend
reports these errors directly.

## Evidence and remaining scope

Tests include exact finite coefficient comparisons with the defining
Pochhammer series, independent HPL Taylor expansion, operator identities,
round-trip serialization and selected numeric secondary checks. Finite tests
are not proof of every high-order/rank performance claim. The J4 endpoint
exchange proof applies to J4, not arbitrary integrals.

Uncovered parameter patterns, nonterminating ranks other than `pF(p-1)` and
arbitrary contour winding remain outside the supported contracts. The
[port audit](port-audit.md) records the shipped HPL special-value and
simplification families now covered by derived identities. The optional
Wolfram bridge supplies live comparison evidence within that documented scope.

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

## HypExp 2 class correspondence

The upstream base-pattern whitelist uses `(p, half-upper-count, half-lower-count)`.
Native kernel selection now covers its listed classes:
`(2,0,1)`, `(2,1,1)`, `(2,2,1)`, `(2,1,0)`,
`(3,0,1)`, `(3,1,0)`, `(3,1,1)`, `(3,2,2)`, `(3,3,2)`,
`(4,1,1)`, `(4,3,3)`, and all-integer `pF(p−1)`.
Integer shifts are implemented with finite-sum, integral, or contiguous
reductions subject to exactness/resource guards. This correspondence is a
coverage map, not exhaustive validation of every parameter choice or parity
with all upstream special cases, continuation rules, or HPL operations.

## Preserving exceptional input parameters

Use `hypra.hypergeometric(...)` (or `sympy.hyper(..., evaluate=False)`) when
validation must inspect every supplied parameter. SymPy's default constructor
may cancel equal upper/lower entries before Hypra receives the expression;
removed fixed-singular parameters cannot be reconstructed. Visible fixed
nonpositive integer lower parameters are rejected even at zero argument or
with a terminating numerator. Identical moving parameters can be canceled
meromorphically in the regulator; original and planned expressions are both
recorded in provenance.

## Principal complex arguments

### Exact complex points inside the disk

`NativeConfig(argument_domain="principal-disk")` permits exact numerical
hypergeometric arguments with provable `0<|z|<1`, including negative real
points. These identities continue the regular origin germ. Zero and the
unit endpoint at one have their own rules. Other unit-circle points and
symbolic complex-domain assertions fall outside this mode.

### Exact points on the principal cut plane

For exact numerical points outside the unit disk, explicitly select
`NativeConfig(argument_domain="principal-cut-plane")`. The supported domain is
C minus the real ray `[1,∞)`, with the endpoint at one handled separately.
Arguments must be provably finite and either nonreal or real below one.
The canonical maps remain on compatible branches: case A maps to the upper
half-plane, case B to the unit disk, and matched-half maps to the right
half-plane off its positive cut. The identity theorem extends the same
origin-germ identities, without exchanging integrals and regulator limits.
This is not acceptance of general symbolic complex-domain assertions or
unspecified cut boundary values.

### Negative real arguments

For `z=-x` with `x>0`, the hypergeometric germ continued from zero has
no singularity on the integration path. The indicial powers and logarithms
in an intermediate HPL representation can still depend on the side used for
each term. Hypra approaches the negative ray from the upper half-plane
throughout the calculation.

Integer kernels have binary HPL words. With `y=x/(1+x)`, the exact map
`z=y/(y-1)` expresses their mixed HPLs at `0<y<1`, including `z=-1`.
The matched-half coordinate is `+i sqrt(x)`. The single-half-lower and
all-half coordinate is `-sqrt(x/(1+x))`, and the single-half-upper
coordinate is `(1-sqrt(1+x))/(1+sqrt(1+x))`. Both negative coordinates
lie in `(-1,0)` and approach it from above. Negating the HPL letters maps
their values into the positive interval, with the upper-side logarithm
constant for any trailing zeros. The fractional map supports mixed binary
HPLs through weight ten, and the negation map supports weight 24.
Contiguous shifts operate on the same
analytic germs before substitution. Rank-one families, including equal
parameter cancellations, use the binomial identity
`1F0(a;;z)=(1-z)^(-a)` directly.
