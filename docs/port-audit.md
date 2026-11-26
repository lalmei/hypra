# HypExp 2 and HPL 2 port audit

This audit targets the public operations and shipped reduction families in
the pinned HypExp 2 archive dated 19-Apr-2011 and standalone HPL 2 archive
dated 03-Nov-2011. [Source hashes and locations](upstream-source.md) identify
the references. The Python implementations use derived recurrences and
identities; upstream source and tables are not included in the distribution.

## Gaps closed

| Area | Result |
|---|---|
| Expansion library | Out-of-template requests fall back to the native dispatcher. Registration order, sorted bases, specialization, failed upgrades and Laurent pole depth have regression coverage. |
| Binomial expansions | Rank-one functions use the binomial kernel for half-integer and shifted bases, including functions reduced to rank one by cancellation. |
| Negative real arguments | The native hypergeometric expansion continues the origin solution along the negative real axis for supported parameter families. |
| Gauss endpoints | Convergence uses the positive real part of the excess. Exact irrational, symbolic-positive and complex bases now reach the Gamma expansion. |
| Logarithmic integrals | Beta/Gamma coefficient extraction replaces the former weight-seven HPL route. Beta moments support log weight 32 and polynomial degree 4096; mixed-log moments admit explicit principal nonreal arguments. |
| HypExp polylogarithm argument relations | Exact pullbacks cover all six cross-ratio arguments for classical and Nielsen functions, with explicit upper/lower cut sides. |
| HypExp named polylogarithm reductions | Exact rational linear algebra reduces all classical/Nielsen functions through weight five to logarithms, classical polylogarithms and canonical `S32(z)`, including all six cross-ratio arguments. |
| HPL known functions | All shipped patterns are covered, including signed Nielsen families and two-letter quartic words. |
| Ordinary endpoint constants | Double-shuffle, duality and height-one identities reduce all convergent MZVs through weight eight. General Euler, repeated-index, Zagier and prefixed families extend through weight 24. |
| Alternating endpoint constants | Colored shuffle/stuffle and distribution reduce every convergent colored sum through weight six, matching the shipped general table ceiling. The retained primitive is `ES((5,1),(-1,-1)) = HPLs6-zeta(6)`. |
| At-i constants | Exact shuffle and argument-map equations reproduce the shipped named-constant basis through weight six, with one, nine and 23 retained primitives at weights four, five and six. Simple classical words reduce through weight 24. |
| Real-cut continuation | Numerical values support all ternary words through weight 24. Direct symbolic inversion supports weight eight and repeated letters through 24. |
| Taylor series | Exact algebraic/transcendental interior centers and explicit nonreal centers are supported. SymPy series use the enlarged real-center domain. |
| Plus/minus notation | Expression transformations, log extraction, Lyndon reduction, integration and series process `HPLPM` atoms. Exact AST and Wolfram transport preserve the at-i primitives, which have independent numerical evaluation. |
| Scalar HPL algebra | Ordinary and plus/minus HPLs of a commutative argument commute in SymPy products. |
| Values at minus one | `hpl_at_minus_one` evaluates finite ternary endpoints and accepts a cut side for branch-dependent words. SymPy evaluates side-independent finite values directly. |
| Wolfram comparison bridge | Nielsen and classical polylogarithms round-trip with distinct heads. The cache fingerprint includes the upstream endpoint reduction tables. |
| Numerical and cache integration | Nielsen outputs evaluate numerically on their supported branches. Reduction maps are immutable; expansion cache fingerprints include every Python module. |
| Distribution | Wheel and source builds include the new modules and exclude upstream archives and local TeX build artifacts. |

The [parity map](parity-map.md) links these operations to their Python APIs.
Each API documents its domain and resource bounds. Expansions and numerical
continuation fail outside their contracts; symbolic reducers preserve exact
unsupported atoms where documented.

## Verification and review

The initial native baseline passed 1,071 tests. Each implementation received
focused validation and an independent review. Checks include exact shuffle
and differential identities, Gamma derivatives, defining logarithmic integrals,
finite nested sums with tail estimates, and high-precision numerical comparisons.
Review found and fixed missing high-weight constant families and a mutable
cache that could corrupt subsequent exact reductions. It also found that
the former at-i evaluator retained unnecessary multiple polylogarithms in
simple closed forms; the generated canonical basis closes that symbolic gap.

The opt-in live suite executes the real pinned packages in Wolfram. It covers
the integer pole/offset paths and half-integer dispatcher classes, all 2,184
ordinary/plus-minus minimal-set reductions through weight six, all 1,093
ordinary at-i words through weight six, six plus/minus at-i cases, and all
1,053 plus/minus at-i reductions of weights four through six by exact symbolic
comparison. It also checks 54 polylogarithm argument/named-function relations,
mixed weight-seven values on both sides of both real cuts, all 485 convergent
colored endpoint constants through weight six, and primitive transport at
weights four, five and six.

The earlier PR #28 audit passed **1,637 native tests** and all **28 live backend tests**
across the integrated and targeted runs. Ruff lint/format checks, Pyright,
the strict documentation build and offline asset checks passed. Wheel and
source-distribution contents were checked, and an isolated wheel installation
passed native expansion, canonical at-i reduction and exact AST smoke checks.

On 2026-11-26, all **1,702 collected tests passed** with the six follow-up
fixes applied. Python 3.12.14 ran the 1,670 native and offline tests.
Python 3.15.2 ran all 32 live tests against the pinned HypExp 2 and HPL 2
packages. The live comparisons included direct and cancellation-reduced
binomials, negative real arguments, scalar HPL products, and both sides of
finite endpoints at minus one. The Wolfram bridge returned the exact Nielsen
and classical polylogarithm coefficients in the same expansion. Editing a
temporary copy of `h6table.mat` forced a cache miss after a verified hit and
reproduced the exact endpoint value.

```sh
uv run --locked pytest -q -m 'not real_backend' --no-cov
HYPRA_RUN_REAL=1 uv run --locked pytest -q -m real_backend --no-cov
uv run --locked --group quality pyright --pythonpath .venv/bin/python
uv run --locked --group quality ruff check .
uv run --locked --group quality ruff format --check src tests examples scripts
make docs
uv build
```

## Remaining distinctions

The port supplies explicit Python operations instead of Mathematica global
rewrite switches or an implicit cut-side preference. Exact comparisons cover
the shipped minimal and at-i bases; they do not imply Mathematica syntax or
global-evaluation compatibility.

Generic higher mixed HPLs are also retained by the shipped known-function
rules. Arbitrary winding paths are outside the shipped real-region
continuation interface. Numerical tail and stability diagnostics do not certify
floating-point roundoff. Singular endpoints and continuation in parameters
remain separate, explicit contracts.

The general-family identities use the published mathematical results in
[Borwein, Bradley, Broadhurst and Lisoněk, Examples 2.2–2.3](https://arxiv.org/abs/math/9910045).
The full audit supports the listed operation coverage and tested equivalences;
it does not assert exhaustive equality of every upstream symbolic expression.
