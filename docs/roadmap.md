# Native port status

## Implemented

- Exact all-order HPL recursion for 2F1(1,bε;1+cε;z).
- Two exact half-integer families, including telescoping odd-polylogarithm sums.
- Whole-expression Gamma/rational Laurent prefactors and conservative depth planning.
- Symbolic HPL words, derivatives, selected known-word reductions, exact JSON and LaTeX.
- Native API/CLI, provenance/cache, independent symbolic/analytic/numerical checks.
- Official source snapshots and a hashed algorithm map for continued porting.

## Next substantive translation work

1. Integer shifts and contiguous/reduction operators for broader 2F1 families.
2. HypExp Z/S nested-sum algebra and operator composition for general pF(p-1).
3. HypExp 2 half-integer class reductions and rationalized-variable recurrences.
4. HPL shuffle/reduction identities and explicit analytic continuation conventions.
5. Extend exact regression cases and external independent oracles where available.

Unsupported families currently raise, rather than passing through an unexpanded
object or claiming full HypExp parity. Higher-order performance work comes after
correctness of each added family. Preserve attribution and clarify upstream
redistribution permissions before publishing direct code translations.

Mathematica is not required and obtaining a Wolfram setup is not the project goal.
The retained external adapter is optional experimental comparison code only.
