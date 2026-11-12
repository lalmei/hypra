# Binary square and fractional HPL maps

Both APIs require `domain="0<z<1"`, an explicit caller assertion, and expanded
binary words of weight at most ten. The bound limits the multilinear expansion
to at most 1024 words. An original `-1` letter is rejected: its pullback would
introduce poles outside the output HPL alphabet.

## Square

`square_hpl(word, z, domain="0<z<1")` returns `H_word(z**2)` in HPLs at `z`.
For `t=z²`, direct differentiation gives

\[
 dt/t=2\,dz/z,\qquad dt/(1-t)=dz/(1-z)-dz/(1+z).
\]

Expand each letter of the iterated integral by multilinearity: `0 -> 2*0`
and `1 -> 1 - (-1)`, where the latter minus denotes a coefficient, not a
compressed index. The resulting alphabet is ternary. There is no extra origin
constant because analytic words vanish at zero, and the logarithmic base
case is exactly `log(z²)=2 log(z)` on the stated positive interval. Shuffle
identities extend the same construction to trailing-zero words.

## Fractional map

`fractional_hpl(word, z, boundary="upper" | "lower", domain="0<z<1")`
returns `H_word(z/(z-1) +/- i0)`. The boundary is mandatory even when a
particular word has no discontinuity. For `t=z/(z-1)`,

\[
 dt/t=dz/z+dz/(1-z),\qquad dt/(1-t)=-dz/(1-z).
\]

For origin-analytic words ending in `1`, expand tensor letters as
`0 -> 0 + 1`, `1 -> -1` (here `-1` means **minus the letter 1**, not the
letter -1). Both sides solve the same differentiated recursion and vanish at
zero. They therefore agree, continued along the negative real ray, where
these analytic-origin binary words have no branch cut.

Trailing-zero words need extra monodromy constants. Shuffle extraction first
expresses them in analytic-origin words and explicit logarithms. Then impose

\[
 \log(t\pm i0)=\log z-\log(1-z)\pm i\pi
              =H_0(z)+H_1(z)\pm i\pi.
\]

The implementation does not silently use the analytic-origin pullback for a
trailing-zero word. Its boundary selection is the side of the **final negative
argument**, not a direction assigned to `z`.

Tests compare square-map exact local series, fractional differential
recursions and origin constants, the independent dilogarithm identity
`Li₂(z/(z-1))=-Li₂(z)-log²(1-z)/2`, and a mixed trailing-zero discontinuity.
No numerical fit or stored transformation table is used.
