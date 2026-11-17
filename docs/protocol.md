# Exact protocol and safety

The default native backend operates entirely in Python/SymPy. Inputs and outputs
use a whitelisted exact JSON AST for reproducibility and CLI/cache interchange.
Integers are strings, fractions retain numerator/denominator, and floats, unknown
heads, invalid symbol names, duplicate JSON keys and malformed structures fail.
No Python eval or executable user expression strings are used.

Native hypergeometric coefficients come from documented mathematical recurrences,
not stored expected values. The complete expression is expanded deeply enough
for every known Laurent pole, then truncated to the requested inclusive order.
A zero truncation and negative requested orders are handled explicitly. Unknown
valuations, unsupported families/domains and non-Laurent coefficients fail closed.
There is no automatic switch to the optional external adapter.

Outputs record exact AST, order, assumptions, algorithm/source/SymPy versions
and cache status. HPL uses expanded letters {-1,0,1}. The distinct `HPLPM`
head encodes minus, zero and plus with integer codes -1,0,1; those codes do
not denote ordinary signed letters. Unsupported alphabets are not silently
reinterpreted. LaTeX output is generated from exact SymPy objects.

The optional Wolfram transport uses fresh processes, private JSON files and a
fixed head whitelist. The bridge requests HypExp's inclusive expansion order
after budgeting all prefactor poles. HypExp can return a normal expression
without an explicit remainder. In that case the bridge relies on the upstream
inclusive order contract. If an explicit series remainder is present, the
bridge verifies its regulator, center and depth before removing it. HPL
output with `plus` and `minus` letters expands exactly into the ordinary
alphabet. Retained at-i `HPLI` primitives instead cross the boundary as
`HPLPM(word,I)`, preserving the canonical constant basis. The same exact
atoms can be sent back through the bridge and evaluated numerically in Python.
