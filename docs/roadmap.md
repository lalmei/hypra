# Native implementation status and roadmap

## Implemented building blocks

- Integer `pF(p-1)` canonical nested sums and contiguous reductions, including
  moving negative-upper reindexing and singular-lower pole budgeting.
- Matched-half colored nested sums and integer-offset reductions.
- Four canonical half-integer rationalized ODE patterns, plus specialized
  elementary/polylogarithmic identities.
- Rational-HPL integration, definite origin constants, shuffle log extraction,
  binary endpoints/MZVs, branch-explicit negation/reflection, table-free Lyndon
  shuffle-basis reduction and exact local series.
- Arbitrary-word numerical HPL evaluation on the real interval and explicit
  principal slit disk; off-cut radial Taylor continuation with propagated
  truncation estimates and separate high-precision stability checks.
- Real-domain SymPy `HPL.evalf()` integration without eager branch changes.
- Optional hard native-backend deadline and exceptional-parameter validation.
- Exact terminating sums and condition-checked Gauss endpoint rewriting.
- Whole-expression Laurent prefactors, exact serialization, LaTeX, provenance,
  cache and native API/CLI.

See [scope](scope.md) for the precise supported families and restrictions.
These are native Python mathematical implementations, not a literal port of
upstream Mathematica routines. No Mathematica setup is needed or planned as a
prerequisite.

## Next work, in priority order

1. Consolidate dispatch/regression coverage across all implemented kernels,
   endpoint preprocessing and complete-expression pole budgets. Include mixed
   families and cancellation-sensitive products, not only isolated functions.
2. Complete exceptional-parameter and branch coverage within the supported
   HypExp class patterns. The exact hypergeometric API has explicit real-interval and
   restricted principal-complex domains; broader numerical HPL domains do not
   automatically widen that contract.
3. Expand independent exact coefficient benchmarks, higher-rank cases and
   complexity/performance measurements. Optimize algebra without weakening
   fail-closed checks or origin-limit verification.
4. Extend the existing alternating endpoint, height-one MZV and symbolic-map
   implementations with independently proved higher-weight reductions. Binary
   weight-three classical reduction and explicit real shuffle regularization
   are now available; neither is full upstream special-value-table parity.
5. Add explicit boundary/winding-path continuation where needed; the numerical
   API currently supports only its stated off-cut radial path. Broader contour
   integration is not inferred. Divergent endpoint regularization is available
   only through explicitly named shuffle-regularization APIs, never as a limit.

Full HypExp 2/HPL parity is incomplete. A future comparison oracle may help,
but the retained external adapter is experimental, unvalidated and optional.
Keep academic attribution and source-provenance records; upstream licence
restrictions still apply to copying or redistributing upstream code/tables.
