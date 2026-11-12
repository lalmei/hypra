# Inversion with explicit boundary values

`invert_hpl(word, z, boundary="upper" | "lower", domain="0<z<1")`
returns an exact expression for `H_word(1/z + i0)` or `H_word(1/z - i0)`.
The boundary is mandatory and refers to the final argument `1/z`. It is not
inferred from a floating-point sign or from an arbitrary complex assumption.
Ternary words through weight six are supported subject to the component
transformations' resource bounds.

## Composition and branch direction

Write `M(u)=(1-u)/(1+u)`, an involution. Then
`M(-M(z))=1/z`. First express `H_word(M(u))` in HPLs at `u` using the exact
Möbius transformation. Continue `u` from the positive interval through the
upper half-plane to `-M(z)+i0`. Since

\[
 \operatorname{Im} M(u)=\frac{-2\operatorname{Im}u}{|1+u|^2},
\]

its image approaches `1/z` from **below**. The negation transformation's
explicit `log(-v+i0)=log(v)+i*pi` therefore gives the lower final boundary.
A second real Möbius transformation puts all remaining HPLs at `z`.

This path avoids the singular points `0,+1,-1`; near the negative segment
`(-1,0)` all other logarithms retain their real determinations. The original
Möbius identity extends along that path by analytic continuation, with the
constants fixed by its original real-interval endpoint construction.

All final HPLs at `z`, alternating Euler sums, MZVs and elementary constants
are real on the declared interval. Schwarz reflection therefore obtains the
upper boundary by reversing explicit `i`. No conjugation of an unspecified
complex variable or constant is assumed.

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
just numerical samples. Arbitrary variable-dependent functions in an
intermediate expression are rejected rather than silently continued.
Unreduced higher-weight endpoint constants remain exact Euler sums.
