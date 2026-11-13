# Complete real classical reduction at weight three

`classical_hpl(..., domain="0<z<1")` now reduces **all 27 ternary words of
weight three**, not merely words using two letters. This completes weight
three, not arbitrary higher-weight classical reduction.

Let `F_+(z)=H_(0,1,1)(z)` and `F_-(z)=H_(0,-1,-1)(z)`, whose real classical
formulas were derived earlier. The first missing Lyndon primitive is

\[
 A(z)=H_{-1,0,1}(z)=\log(1+z)Li_2(z)
       +\tfrac12F_+(z^2)-F_-(z)-F_+(z).
\]

Integration by parts gives the first product plus
`∫₀ᶻ log(1+t)log(1-t) dt/t`. Use
`2 log(1+t)log(1-t)=log²(1-t²)-log²(1+t)-log²(1-t)` and substitute `u=t²`
in its first term. This proves the displayed formula with zero origin
constant. All arguments stay on the real principal branches.

For the second primitive `B(z)=H_(-1,1,0)(z)`, let
`M(z)=(1-z)/(1+z)`, `h=H_-1(z)`, and `g=H_(-1,1)(z)`. The already-proved
Möbius pullback, followed by exact shuffle reduction, gives

\[
 A(M(z))=C-\zeta(2)h-\tfrac12\log2\,h^2+\tfrac16h^3
         -\log2\,g+gh-H_{-1,-1,0}(z)-2H_{-1,-1,1}(z)-B(z).
\]

Everything except `B` is already classical. Solving this linear equation
therefore determines `B` without a new integration table. The constant is
the convergent endpoint of the preceding `A` formula:

\[
 C=A(1)=\log2\,\zeta(2)-\tfrac32\zeta(3)
      +Li_3(1/2)+\log2\,Li_2(1/2)+\tfrac13\log^3 2.
\]

Retaining `Li₃(1/2)` avoids relying on an additional special-value identity.
`M(z)` and every argument inside the substituted real formulas remain within
`(0,1)` where required. Thus there are no hidden complex-log cancellations.

Lyndon shuffle reduction expresses each remaining ternary word in these two
primitives and the already-supported two-letter words. Tests check the exact
square-log differential construction and zero origin, the exact Möbius/shuffle
equation fixing `B`, and all six genuinely ternary words against independent
HPL numerical evaluation at `1/3` and `3/4` to 29 digits. Previous tests cover
the other 21 words. Weight-four mixed words remain lossless HPLs.
