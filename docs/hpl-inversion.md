# Inversion with explicit boundary values

`invert_hpl(word, z, boundary="upper" | "lower", domain="0<z<1")`
returns an exact expression for `H_word(1/z + i0)` or `H_word(1/z - i0)`.
The boundary is mandatory and refers to the final argument `1/z`. It is not
inferred from a floating-point sign or from an arbitrary complex assumption.
Ternary words through weight eight are supported. Repeated-letter words
extend to weight 24.

## Pulled-back differential equations

Under `x=1/z`, the three kernel forms transform as

\[
 f_0(x)dx=-f_0(z)dz,\qquad
 f_1(x)dx=(f_0(z)+f_1(z))dz,\qquad
 f_{-1}(x)dx=(-f_0(z)+f_{-1}(z))dz.
\]

Starting with the empty word, the implementation integrates the transformed
suffix recursively by prepending these letters. At `z→1-`, the endpoint
logarithm of `H_word(1/z±i0)` is `T±i*pi+log(z)`, where
`T=-log(1-z)`. The last term vanishes. The exact endpoint polynomial fixes
the integration constant by matching its constant term after this shift.
This also handles divergent leading-one words without assigning them finite
ordinary endpoint values. Repeated letters use powers of the weight-one
formula divided by the word factorial.

## Independent low-weight checks

For `s=+1` on the upper boundary and `s=-1` on the lower boundary,

\[
\begin{aligned}
 H_0(1/z+s i0)&=-\log z,\\
 H_1(1/z+s i0)&=\log z-\log(1-z)+s i\pi,\\
 H_{-1}(1/z+s i0)&=\log(1+z)-\log z,\\
 Li_2(1/z+s i0)&=-Li_2(z)+\pi^2/3-\tfrac12\log^2z-s i\pi\log z.
\end{aligned}
\]

These identities test both the integration constants and the branch signs.
Additional mixed-word tests check the exact differentiated chain rule, not
just numerical samples. Classical inversions through weight eight are also
checked against the Bernoulli-polynomial inversion identity.
Unreduced higher-weight endpoint constants remain exact Euler sums.
