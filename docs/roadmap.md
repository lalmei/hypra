# Native implementation status and roadmap

## Implemented building blocks

- Integer `pF(p-1)` canonical nested sums and contiguous reductions, including
  moving negative-upper reindexing and singular-lower pole budgeting.
- Matched-half colored nested sums and integer-offset reductions.
- Three canonical half-integer rationalized ODE systems, plus specialized
  elementary/polylogarithmic identities.
- Rational-HPL integration, definite origin constants, shuffle algebra and
  independent local-series checks.
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
2. Extend half-integer class reductions beyond matched-half patterns and the
   three canonical ODE systems. Derive missing shifts with explicit boundary
   conditions and branch choices; reject unsupported patterns in the meantime.
3. Expand independent exact coefficient benchmarks, higher-rank cases and
   complexity/performance measurements. Optimize algebra without weakening
   fail-closed checks or origin-limit verification.
4. Add justified endpoint/special-value identities beyond Gauss and termination.
   Do not reuse a meromorphic continuation as a divergent endpoint value.
5. Extend HPL reductions and explicit continuation domains/paths. General
   contour continuation and generic integration remain separate projects.

Full HypExp 2/HPL parity is incomplete. A future comparison oracle may help,
but the retained external adapter is experimental, unvalidated and optional.
Keep academic attribution and source-provenance records; upstream licence
restrictions still apply to copying or redistributing upstream code/tables.
