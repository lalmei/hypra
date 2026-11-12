# Real Möbius transformation of HPL arguments

`hypra.hpl_mobius.mobius_hpl(word, z, domain="0<z<1")` expresses
`H_word((1-z)/(1+z))` in HPLs at `z`. Every ternary word through weight six
is supported, subject to the integration engine's resource limits. The domain
is an explicit caller assertion. The implementation makes no claim about
complex continuation or inversion across a branch cut.

Put `t=(1-z)/(1+z)`. Both `t` and `z` lie in `(0,1)`. Direct differentiation
of the three defining one-forms gives

\[
 \frac{dt}{t}=-\frac{dz}{1-z}-\frac{dz}{1+z},\qquad
 \frac{dt}{1-t}=-\frac{dz}{z}+\frac{dz}{1+z},\qquad
 \frac{dt}{1+t}=-\frac{dz}{1+z}.
\]

Thus the alphabet closes without introducing new singular points. For a word
with first letter `0` or `-1`, its endpoint at `t=1` is convergent; the
transformed kernel is regular at `z=0`. The recursion uses that exact endpoint
plus the integral from zero to `z` of the pulled-back kernel times the already
transformed tail. Tails have at most logarithmic singularities at zero, so
these integrals converge. Rational-HPL integration evaluates them exactly.

Endpoint constants come from `hpl_alternating_at_one`, not an identity table.
Unreduced constants remain exact `EulerSum` objects with strict decreasing
indices, explicit colors and the endpoint evaluator's outer-cutoff convention.
Consequently the transform does not assume an unproved reduction to powers of
pi, ordinary zeta values or logarithms.

A leading `1` can make the endpoint divergent. Instead of assigning a
regularized value, shuffle `H_1 H_tail` to isolate the word. Other same-weight
terms have fewer leading ones, so this recursion terminates. The base case is

\[
 H_1(t)=-\log2-\log z+\log(1+z),\qquad
 H_{1,\ldots,1}(t)=H_1(t)^n/n!.
\]

There is no branch ambiguity in these real logarithms. The transformation is
an involution at the function level. Simplifying a symbolic double transform
may require additional Euler-sum identities; a syntactically nonzero constant
remainder is not itself evidence of a failed transformation.

Validation checks differential pullbacks for all 39 nonempty ternary words
through weight three and selected weight-four/five words. Independent
elementary weight-one constants, repeated-letter shuffle identities, domain
requirements and resource limits are also checked. Constant determination is
backed by the separate alternating-endpoint evaluator's tests.
