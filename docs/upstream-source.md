# Upstream source study: HypExp 2 and HPL 2

## Acquired references, not runtime dependencies

On 12 November 2026, the official University of Zurich distribution archives were
downloaded into local `reference/upstream/` and extracted with path/link checks.
**No upstream installer or Mathematica source was executed.** Extracted files are
read-only. They are excluded from Git and from the Python distribution.
The exact URLs, byte counts and SHA-256 hashes of every downloaded source file
are in [the source manifest](upstream-source-manifest.json).

| Reference | Version marker in archive | Archive SHA-256 |
|---|---|---|
| [HypExp 2.0](https://www.physik.uzh.ch/data/HypExp/HypExp-2.0.tar.gz) | 19-APR-2011 | `d6275525a37d393c7706210d875364d06962f84e1715b570f07c56ed3cf9a60f` |
| [HPL 2.0](https://www.physik.uzh.ch/data/HPL/HPL-2.0.tar.gz) | 03-NOV-2011 | `ca965edb320c735f31433daf11b099ef11d2fbd1a8108b11293d18a6674136ee` |

The HypExp archive bundles an older HPL revision dated 28-SEP-2010. It is **not**
byte-identical to the separately downloaded HPL revision. Keep this distinction
in any future reference comparison.

## Licence and attribution

The [HypExp repository record](https://elsevier.digitalcommonsdata.com/datasets/cphpktd5cn/1)
and [HPL 2 repository record](https://elsevier.digitalcommonsdata.com/datasets/hntx3r95cx/1)
identify the CPC licence. The [official licence](https://www.elsevier.com/about/policies-and-standards/open-access-licenses/elsevier-user/cpc)
permits academic/non-profit use within the specified research-group scope,
requires the relevant CPC publication to be cited, and requires author permission
for third-party redistribution of original or modified program code. It is not
an MIT/BSD-style permission to redistribute a port. A file called `HPLnolicense`
in the archives is an installation placeholder, **not** a public-domain grant.

Accordingly, this project keeps the acquired source local for study, does not
copy upstream tables into the wheel, and uses independently written Python
implementations of documented mathematical identities. Before distributing a
literal translation or derivative of upstream code, obtain the authors' written
permission and settle the applicable licence. Academic purpose alone does not
remove the redistribution restriction. This is a provenance safeguard, not a
legal opinion.

## Concrete source map for translation work

The following line numbers refer to the pinned archives in the manifest; they
are an inspection map, not a claim that these routines have already been ported.

| Component | Source location | Python translation target |
|---|---|---|
| Parameter classification | `HypExp.m`, lines 1879–1942 | Exact affine parameter parsing, integer/half-integer classification, reject unsupported classes explicitly. |
| Gamma/Pochhammer expansion | `HypExp.m`, `GammaExpand` near line 300 | Rational truncated coefficient algebra, exact gamma prefactors, and pole-aware order planning. |
| Integer-parameter nested-sum engine | `HypExp.m`, `ZPart` near 1708 and `MyHypExpansion2` 1721–1806 | Expand coefficient ratios into finite nested sums; reduce shifted denominators; map sums to HPLs with exact coefficients. |
| Half-integer reduction | `HypExp.m`, `HaHypExpAlg` 1567–1588 and `HaHypExp` 1626 onward | Reduce parameter offsets using derivative/integral operators and explicitly supported base families; track transformed arguments/branches. |
| Half-integer family dispatch | `HypExp.m`, `JrsKnown` line 1623 | Explicit supported-family registry rather than an unrestricted claim about every half-integer pFq. |
| Unit-argument path | `HypExp.m`, 2004–2023 and `HypExpAtOne` near 2418 | Separate endpoint treatment; do not substitute z=1 into a generic expansion without convergence checks. |
| HPL index notation | `HPL.m`, 163–250 | Lossless compressed/expanded index conversion with explicit `0, 1, -1` alphabet and trailing-zero handling. |
| HPL product algebra | `HPL.m`, 394–432 | Shuffle recursion with multiplicities; avoid factorial permutation enumeration where dynamic programming suffices. |
| HPL argument transformations | `HPL.m`, 1452 onward | Branch-aware transformations and boundary constants; keep unsupported continuation explicit. |

### Historical implementation sequence

The initial sequence below has largely been implemented; it is retained as
design provenance, not a current list of unimplemented tasks. The
[parity map](parity-map.md) is the current capability audit.

1. Implement a narrow integer-parameter family using its differential equation
   and exact iterated-integral recursion. Validate coefficients independently
   against the defining power series, not against hard-coded acceptance outputs.
2. Add the HPL shuffle/derivative primitives needed by that recursion, retaining
   symbolic HPLs whenever a reduction is not proved.
3. Extend contiguous shifts and generic nested-sum denominators incrementally.
4. Add a documented half-integer base family and its branch map before broader
   half-integer dispatch. Do not equate a cancelled elementary half-integer case
   with general HypExp 2 parity.
5. Treat z=1 limits, endpoint integration, MZV reduction and numerical continuation
   as distinct features requiring their own convergence and branch tests.

## HPL 2.0 functionality audit

This audit uses the **standalone 03-Nov-2011 HPL.m**, the manifest-pinned
archive, and the public abstract of
[Extension of HPL to complex arguments](https://arxiv.org/abs/hep-ph/0703052).
The bundled `HPL_README` ends with a September 2005 date, so it is not by
itself a reliable description of the 2.0 additions. No source was executed,
no table entries were copied and no installer was run.

| Evidence in pinned HPL.m | What it establishes | Native status / distinction |
|---|---|---|
| Public declarations, lines 31–101 | Product/log extraction, known-function/minimal-set reduction, argument changes, integration, index conversions, complex continuation, plus/minus notation, and HPLI | Existing Python APIs cover these operation categories with narrower explicit contracts; a shared category is not full parity |
| Continuation region and side declarations, lines 61–69 | Real exterior/negative intervals and a controlled approach side; upstream also has a default side | Native symbolic region maps and numerical boundary values require an explicit side; arbitrary winding paths remain unsupported |
| Integration declaration, lines 73–77 | Rational kernels with poles at 0 and ±1 multiplying HPLs | Native rational-HPL primitives and origin-normalized integration check unsupported poles and divergent lower limits |
| Argument-transform declaration, lines 79–86 | Negation, reflection, square-root, fractional, Möbius and inversion transformations, with branch-side data | Native maps are separately proved and resource-bounded; no generic complex automatic rewrite is implied |
| Plus/minus and HPLI declarations, lines 89–101 | Extended kernel notation and values at the imaginary unit | Native multilinear notation conversion; exact at-i low-weight/simple-word reductions, not a full table import |
| Optional minimal-set loaders, lines 685–697 | Basis rules depend on available weight-specific files | Native Lyndon reduction uses an algorithmic basis; upstream table-dependent basis equivalence is unverified |
| Optional MZV/at-i table loaders, around 1804 and 1828–1833 | Higher special-value coverage depends on separately loaded data | Native exact remainder constants and bounded numerics do not reproduce those named tables |

The at-i work derives all weight-two ternary values from iterated integrals
and shuffle identities, retaining classical dilogarithms when appropriate.
Numerical real-cut evaluation composes existing exact region maps with interior
numerics, rather than approximating the cut by an arbitrary imaginary offset.
These are native mathematical implementations, not execution of HPL 2.0.

The next bounded symbolic gap is the remaining mixed-sign weight-three
classical primitives, followed by broader constant relations and justified
continuation domains. See [the current parity map](parity-map.md) for the
specific proposed basis and the separation of implemented work from plans.

### Distribution boundary remains unchanged

The original HPL 2 archive, Mathematica routines, minimal-set files, at-i tables
and numerical constant tables remain local ignored reference material. The CPC
licence record cited above does not become permissive merely because the goal
is an academic Python package. Independently derived algorithms need their own
implementation and attribution; any direct source/table translation or
redistribution requires a separate permissions review. No claim of legal
clearance for an arbitrary literal port follows from this technical audit.

## Important observation about the abandoned wrapper-first approach

The inspected public HypExp routines call `Normal[Series[...]]` in several paths
(e.g. lines 1945, 1985, 2047 and 2056), discarding Mathematica's explicit O-term.
A wrapper that insists on seeing `SeriesData` cannot validate those outputs and
must reject them. Source inspection does **not** retroactively turn mock transport
tests into a successful HypExp execution. The current user-directed goal is a
native Python implementation, not installing Mathematica to work around this.

## Citations

- T. Huber and D. Maître, *HypExp, a Mathematica package for expanding
  hypergeometric functions around integer-valued parameters*, Comput. Phys.
  Commun. 175 (2006) 122–144, [hep-ph/0507094](https://arxiv.org/abs/hep-ph/0507094),
  [DOI 10.1016/j.cpc.2006.01.007](https://doi.org/10.1016/j.cpc.2006.01.007).
- T. Huber and D. Maître, *HypExp 2, Expanding Hypergeometric Functions about
  Half-Integer Parameters*, Comput. Phys. Commun. 178 (2008) 755–776,
  [arXiv:0708.2443](https://arxiv.org/abs/0708.2443),
  [DOI 10.1016/j.cpc.2007.12.008](https://doi.org/10.1016/j.cpc.2007.12.008).
- D. Maître, *HPL, a Mathematica implementation of the harmonic polylogarithms*,
  Comput. Phys. Commun. 174 (2006) 222–240,
  [hep-ph/0507152](https://arxiv.org/abs/hep-ph/0507152),
  [DOI 10.1016/j.cpc.2005.10.008](https://doi.org/10.1016/j.cpc.2005.10.008).
- D. Maître, *Extension of HPL to complex arguments*, Comput. Phys. Commun. 183
  (2012) 846–853, [hep-ph/0703052](https://arxiv.org/abs/hep-ph/0703052),
  [DOI 10.1016/j.cpc.2011.11.015](https://doi.org/10.1016/j.cpc.2011.11.015).
