# Branch-explicit HPL argument transformations

`hypra.hpl_transform` provides `negate_hpl(word, z, domain="0<z<1")`
and `reflect_hpl(word, z, domain="0<z<1")`. Each returns an exact SymPy
expression for the indicated HPL, not a numerical approximation. The domain
argument is a caller assertion: the implementation does not infer it from an
unconstrained symbol. No transformation is applied automatically to a backend
result with unspecified branches.

## Negation

For a word with nonzero final letter, substitution of the integration path
`t=-u` gives

\[
 H_{a_1,\ldots,a_n}(-z)=(-1)^{N_{\ne0}}
 H_{-a_1,\ldots,-a_n}(z).
\]

Indeed `dt/t=du/u`, whereas both nonzero kernels change sign and exchange
letters `1` and `-1`. The nonzero final letter ensures an analytic origin
and fixes the integration constant to zero. This formula alone is **not**
valid for a trailing-zero word.

For trailing zeros, the shared shuffle-log extraction recursively expresses
the word in analytic-origin words and powers of `H_0`. The principal value
on the **upper boundary** of the negative-real cut is explicitly
`H_0(-z)=log(z)+i*pi`. Applying that identity to the extracted logarithms
includes the required monodromy. A lower-cut approach would instead require
`-i*pi` and is not this API's convention. Weight is capped at 24, consistent
with the underlying shuffle engine.

## Reflection

Reflection currently supports every word over `{0,1}` through weight eight.
This is a resource bound, not a restriction to a short table of known
identities. Letter `-1` would generate a pole at `z=2`, outside the HPL
alphabet, and is rejected.

For a word beginning with zero, differentiation and its convergent value at
one give

\[
 H_{0,w}(1-z)=H_{0,w}(1)
 -\int_0^z\frac{H_w(1-t)}{1-t}\,dt.
\]

The tail is recursively reflected; exact rational-HPL integration evaluates
the integral. Reflected tails have at most logarithmic singularities at zero,
so this integral converges. `hpl_at_one` supplies the endpoint from its
geometric-series/MZV construction and shuffle reduction; there is no copied
endpoint table here. Unreduced constants remain exact `MZV` expressions.

A word beginning with one is first reduced by shuffling `H_1 H_tail`.
Other same-weight terms have strictly fewer leading ones, so the recursion
terminates. Pure-one words use
`H_(1,...,1)(1-z)=(-log(z))^n/n!`. These steps also cover words whose endpoint
at one diverges: no divergent endpoint is assigned a finite value.

The implementation fixes constants using actual endpoint values, not just
matching derivatives. Tests check all binary words through weight three
by differential equations and double reflection, exact origin series for
negation, explicit trailing-zero monodromy, selected weight-four/five
reflection derivatives, and rejection of unsupported domains/letters/weights.
It is not a general minimal-basis reducer or a continuation-path engine.

## Whole-expression argument unification

`transform_hpl_arguments(expr,z,domain="0<z<1",boundary="upper")` dispatches
recognized HPL arguments to a common z: −z, 1−z, z², (1−z)/(1+z), z/(z−1),
and 1/z. It requires an explicit boundary for cut-valued arguments and rejects
unknown variable-dependent arguments. Constant-argument HPLs remain unchanged.
The component transforms' alphabet and weight restrictions still apply; this
wrapper does not infer new branch identities from syntactic resemblance.
