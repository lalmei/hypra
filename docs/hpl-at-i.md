# Exact HPL values at i

`hpl_at_i(word)` returns the continuation from the origin through the upper
quadrant to `i`, with the logarithmic base case `log(i)=i*pi/2`. It supports
all ternary words through weight six, plus repeated single letters and
classical words `0^(n-1),±1` through weight 24. Mixed words use named constants
and the retained `HPLPM(word, I)` primitives of the shipped HPL 2 basis. Use
`hpl_pm_at_i(word)` for plus/minus words through weight six.

## Branch justification

The generator identities continue from `(0,1)` into the upper quadrant.
The independent low-weight evaluator also uses the proved integral identities
from `classical_hpl`, continued along `z=i*t`, `0<t<=1`.
The logarithms `log(z)` and `log(1±z)` keep their principal values along this
path. The dilogarithm arguments `±z` and `(1±z)/2` never cross `[1,infinity)`.
Consequently the continued formulas preserve their constants and branches.
This continuation does not change the real-only public contract of
`classical_hpl`.

In particular,

\[
 H_0(i)=i\pi/2,\quad
 H_1(i)=-\tfrac12\log2+i\pi/4,\quad
 H_{-1}(i)=\tfrac12\log2+i\pi/4.
\]

For the dilogarithm, splitting its convergent series into even and odd powers
gives

\[
 Li_2(\pm i)=-\pi^2/48\pm iG,
\]

where `G` is Catalan's constant. The real even sum is `-ζ(2)/8`; the imaginary
odd sum is the defining Dirichlet beta value `β(2)`. The reflection identity
at `1/2` gives `Li₂(1/2)=π²/12-log²2/2`. These are derived reductions, not
numerical fits.

Landen's dilogarithm identity reduces the mixed weight-two constants to
`pi`, `log(2)` and Catalan's constant. All 27 ternary cubic words use those
constants, `zeta(3)`, and `Li3((1-i)/2)`. The square-log identity used for
cubic words reaches `1-z²` on its lower polylogarithm cut along `z=i*t`;
the implementation keeps that side when simplifying the constants.

Repeated-letter words use powers divided by factorials. Classical higher
words at `±i` split into their even Dirichlet-eta and odd Dirichlet-beta
series. Odd beta orders reduce to Euler numbers and powers of pi; even
orders use polygamma at `1/4`.

Tests independently integrate the defining one-dimensional path integral for
all convergent weight-two cases to 35 digits, using the explicit logarithmic
base case for `H_00`. Exact principal-log/Catalan identities, higher simple
words, and unsupported-input rejection are also checked. The numerical
comparisons supplement the branch and integral derivation.

## Mixed words at weights four through six

The canonical reducer constructs a finite rational linear system at each
weight. Its equations are shuffle products, the square substitution `u=z²`,
the Möbius substitution `M(z)=(1-z)/(1+z)`, and the classical values at `±i`.
Fractional-argument pullbacks supply `Li3((1-i)/2)` and `Li4((1-i)/2)`.
Endpoint logarithmic polynomials fix the integration constants.

Continuing the Möbius identity through the upper quadrant reaches `M(i)=-i`
through the lower quadrant. Analytic negation changes the sign of each plus
letter; trailing zeros use the explicit logarithmic shift `-i*pi`. Under the
square substitution, zero kernels acquire a factor of two and the minus
kernel becomes the ordinary one kernel. The endpoint is `-1` with
`log(-1)=i*pi` on the induced upper boundary.

These equations leave exactly one primitive at weight four, nine at weight
five, and 23 at weight six. Their word names select the coordinates used by
the pinned HPL tables; all coefficients are derived by elimination. No
formula table or numerical relation fitting is used. Both ordinary and
plus/minus input words reduce to the same basis. For example,

\[
 H_{+,0,0,+}(i)=2G^2-\frac{\pi^4}{32}.
\]

The retained constants are `HPLPM(word, I)` objects. They round-trip through
the exact JSON AST and have independent numerical evaluation, including
when returned by the optional Wolfram bridge.

### Independent numerical evaluation

For a retained primitive, split the straight path at `i/2`. The path composition formula sums over
all ways of dividing the word between the two segments, with the prefix
on the later segment. Evaluate that segment in reverse from `i` to
`i/2` and reverse the prefix, multiplying by `(-1)^prefix_length`.
See [Duhr, Gangl and Rhodes](https://arxiv.org/abs/1110.0458) for the
iterated-integral and multiple-polylogarithm conventions.

For a segment with displacement `d`, subtract its start from the kernel
poles. For the nonzero translated poles `a1,...,ak`, the multiple
polylogarithm arguments are `d/a1,a1/a2,...,a(k-1)/ak`. Zero poles give
the usual leading-zero blocks in the weight indices. On the first
segment, shuffle extraction removes trailing zeros and fixes the
tangential logarithm as `log(i/2)`. In both segment expansions each
cumulative product has modulus at most `1/2`. Geometric kernel
expansion and termwise integration therefore give absolutely convergent
strict nested sums. This numerical path bypasses the canonical reducer, so
evaluating a retained primitive does not recursively ask for its own reduction.

The live Wolfram test compares all 1,093 ordinary words of weights zero through six
and selected plus/minus words against HPL 2.0 at 35 decimal digits. A separate
exact comparison checks all 1,053 plus/minus reductions of weights four through
six against the shipped formulas, after normalizing the conjugate trilogarithm
pair and quarter-period Hurwitz zeta/polygamma expressions into the same
named-constant basis.
