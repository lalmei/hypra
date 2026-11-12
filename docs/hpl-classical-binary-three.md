# Classical reduction of binary weight-three HPLs

`classical_hpl(expression, z, domain="0<z<1")` now reduces every binary
weight-three word to real logarithms and classical polylogarithms. Existing
weight-two ternary reductions and simple all-weight words remain supported.
General higher-weight or mixed ternary words are retained losslessly.

No weight-three table is needed. Shuffle extraction removes leading `1`s and
trailing `0`s. At weight three the only remaining binary words are `001` and
`011`. The first is `Li₃(z)` by its defining iterated integral. For the second,

\[
 H_{0,1,1}(z)=\frac12\log z\,\log^2(1-z)
 +\log(1-z)Li_2(1-z)-Li_3(1-z)+\zeta(3).
\]

Differentiating the right-hand side cancels all terms except
`log²(1-z)/(2z)`, exactly the defining derivative `H_11(z)/z`.
As `z→0+`, both logarithmic products vanish and `Li₃(1-z)→ζ(3)`, fixing
the integration constant to zero. This proves the identity on `(0,1)`.
Applying the already-verified shuffle identities then proves the other six
binary words, rather than introducing separately stored formulas.

Tests check the exact differentiated recursion for all eight words, the new
formula's origin constant, and independent HPL numerical evaluations at
`z=1/3` and `z=3/4` to 29 digits. Numerical checks are secondary evidence;
the derivative, endpoint, and shuffle arguments establish the identities.
