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
There is no automatic switch to the experimental external adapter.

Outputs record exact AST, order, assumptions, algorithm/source/SymPy versions
and cache status. HPL uses expanded letters {-1,0,1}; unsupported alphabets are
not silently reinterpreted. LaTeX output is generated from exact SymPy objects.

The optional Wolfram transport uses fresh processes, private JSON files and a
fixed head whitelist. The bridge requests HypExp's inclusive expansion order
after budgeting all prefactor poles. HypExp can return a normal expression
without an explicit remainder. In that case the bridge relies on the upstream
inclusive order contract. If an explicit series remainder is present, the
bridge verifies its regulator, center and depth before removing it. HPL
output with `plus` and `minus` letters expands exactly into the ordinary
alphabet before it crosses the JSON boundary.
