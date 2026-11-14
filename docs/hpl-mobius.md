# Real Möbius transformation of HPL arguments

`hypra.hpl.mobius.mobius_hpl(word, z, domain="0<z<1")` expresses
`H_word((1-z)/(1+z))` in HPLs at `z`. Every ternary word through weight eight
is supported. The domain
is an explicit caller assertion. The implementation makes no claim about
complex continuation or inversion across a branch cut.

Put `t=(1-z)/(1+z)`. Both `t` and `z` lie in `(0,1)`. Direct differentiation
of the three defining one-forms gives

\[
 \frac{dt}{t}=-\frac{dz}{1-z}-\frac{dz}{1+z},\qquad
 \frac{dt}{1-t}=-\frac{dz}{z}+\frac{dz}{1+z},\qquad
 \frac{dt}{1+t}=-\frac{dz}{1+z}.
\]

The transformed letters stay in the same alphabet. Prepending each
pulled-back letter to the transformed suffix gives a primitive with the
required derivative. The integration constant is fixed by the exact
endpoint polynomial of the original word. If `T=-log(1-t)`, then
`T=-log(z)-log(2)+O(z)` near the transformed origin. The constant term
therefore equals the endpoint polynomial evaluated at `T=-log(2)`.
Every nonempty HPL in the primitive has zero constant term in its
origin power/log expansion.

Unreduced endpoint coefficients remain exact `EulerSum` or `MZV`
objects. A divergent HPL is matched by its full logarithmic asymptote.
The transformation is
an involution at the function level. Simplifying a symbolic double transform
may require additional Euler-sum identities; a syntactically nonzero constant
remainder is not itself evidence of a failed transformation.

Validation checks differential pullbacks for all 39 nonempty ternary words
through weight three and selected weight-four/five words. Independent
elementary weight-one constants, repeated-letter shuffle identities, domain
requirements and resource limits are also checked. Constant determination is
backed by the separate alternating-endpoint evaluator's tests.
