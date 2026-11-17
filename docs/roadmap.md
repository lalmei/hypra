# Native implementation status and roadmap

## Implemented building blocks

- Integer `pF(p-1)` canonical nested sums and contiguous reductions, including
  moving negative-upper reindexing and singular-lower pole budgeting.
- Matched-half colored nested sums and integer-offset reductions.
- Four canonical half-integer rationalized ODE patterns, plus specialized
  elementary/polylogarithmic identities.
- Rational-HPL integration, definite origin constants, shuffle log extraction,
  ordinary MZVs through weight eight and colored sums through weight six,
  branch-explicit argument maps, table-free Lyndon reduction and exact local
  series at real or explicitly selected complex centers.
- Classical and Nielsen argument transformations at all six cross ratios,
  named-function reductions through weight five, and the shipped HPL
  known-function families, including signed Nielsen and quartic patterns.
- Arbitrary-word numerical HPL evaluation through weight 24 on the real
  interval, explicit principal slit disk and both sides of the real cuts;
  off-cut Taylor continuation with propagated truncation estimates and
  separate high-precision stability checks.
- Principal off-cut complex and interior real SymPy `HPL.evalf()`.
- Optional hard native-backend deadline and exceptional-parameter validation.
- Exact terminating sums, convergent Gauss and supported higher-rank unit
  endpoints, and opt-in meromorphic continuation in the parameters at `z=1`.
- Exact beta-log moments through log weight 32 and polynomial degree 4096;
  mixed-log definite integrals with real or explicit principal nonreal arguments.
- Whole-expression Laurent prefactors, exact serialization, LaTeX, provenance,
  cache and native API/CLI.

See [scope](scope.md) for the precise supported families and restrictions.
These are native Python mathematical implementations, not a literal port of
upstream Mathematica routines. No Mathematica setup is needed or planned as a
prerequisite.

## Port verification and further work

The [port audit](port-audit.md) records the source inventory, closed gaps,
independent mathematical checks and live Wolfram comparisons. The optional
external adapter has been exercised against the pinned packages. It remains
a comparison tool; native computation does not require Wolfram.

Mixed at-i values of weights four through six now use the shipped primitive
basis, derived by exact shuffle and argument-map equations. Their symbolic
comparison is separate from the independent numerical checks recorded in
the audit.

Further extensions can raise explicit resource ceilings, improve performance
or support arbitrary winding paths. These exceed the bounded operation
contracts in the current [parity map](parity-map.md). New reductions should
continue to use proved identities, exact regression checks and independent
numerical evidence.

Keep academic attribution and source-provenance records. Upstream licence
restrictions apply to copying or redistributing upstream code and tables.
