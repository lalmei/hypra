# Additional mixed-sign weight-three classical reductions

**Update:** All genuinely ternary weight-three words are now also reduced; see [the complete construction](hpl-classical-ternary-three.md). The scope below records this intermediate derivation.

`classical_hpl(..., domain="0<z<1")` now reduces every weight-three word
using at most two distinct letters of `{-1,0,1}`. This adds all words over
`{-1,1}`. The two remaining genuinely ternary Lyndon primitives
`(-1,0,1)` and `(-1,1,0)` are not claimed reduced; words involving all three
letters remain HPLs.

## Real primitives

Define

\[
 K(u)=-\tfrac12(\log2+\log u)^2\log(1-u)
      -(\log2+\log u)Li_2(u)+Li_3(u).
\]

Direct differentiation gives
`K'(u)=(log2+log u)²/[2(1-u)]`. On the declared interval, both
`u_-=(1-z)/2` and `u_+=(1+z)/2` lie strictly between zero and one.
Substitution in the defining HPL integrals yields

\[
 H_{-1,1,1}(z)=Li_3(1/2)-K(u_-),\qquad
 H_{1,-1,-1}(z)=K(u_+)-Li_3(1/2).
\]

The **opposite signs** follow from the opposite directions of the two
substitutions. `K(1/2)=Li₃(1/2)` fixes the origin constants. Every logarithm
and polylogarithm remains real; no boundary continuation is used.

The other mixed-sign Lyndon primitive follows from shuffle alone:

\[
 H_{-1,-1,1}=
 \tfrac12 H_{-1}^2H_1-H_{-1}H_{1,-1}+H_{1,-1,-1}.
\]

Together with `H_(-1,1,1)`, these give all two-letter words by abstract
shuffle extraction. The derivation introduces no copied table or numerical
fit. Euler's dilogarithm reflection identity at `u_±` can simplify equivalent
expressions arising from different shuffle orders.

Tests check exact primitive derivatives and origin constants, every two-letter
word against independent HPL evaluation at `1/3` and `3/4` to 29 digits, and
preservation of the unreduced genuinely ternary primitives.
