# Exact HPL logarithm extraction

`extract_hpl_logs(expr)` removes leading 1 letters and trailing 0 letters by
shuffle identities, retaining explicit powers of `-log(1-z)` and `log(z)`.
For example H(1,0;z) = -log(z)log(1-z) - Li₂(z).

This preserves the argument and does not choose a new continuation path.
Repeated identical letters reduce to the corresponding weight-one HPL power
divided by the factorial. General remaining words stay exact; no claim of a
minimal transcendental basis is made. Each shuffle retains multiplicity.

The recursive trailing-zero reduction decreases the number of trailing zeros
at fixed weight; leading-one removal similarly decreases leading ones. Lower
weight factors terminate the recursion. Identities are independently tested
using local exact HPL series, including mixed-sign words.
