# Exact HPL values at i

`hpl_at_i(word)` returns the continuation from the origin through the upper
quadrant to `i`, with the logarithmic base case `log(i)=i*pi/2`. It supports
all ternary words through weight four, plus repeated single letters and
classical words `0^(n-1),±1` through weight 24. Mixed cubic and quartic
values can include exact multiple polylogarithm constants.

## Branch justification

The weight-two identities are the same proved integral identities used by
`classical_hpl` on `(0,1)`. Continue them along `z=i*t`, `0<t<=1`.
The logarithms `log(z)` and `log(1±z)` keep their principal values along this
path. The dilogarithm arguments `±z` and `(1±z)/2` never cross `[1,infinity)`.
Consequently the continued formulas preserve their constants and branches.
The implementation first obtains the real formula and then makes this
justified continuation; it does not claim the real-only classical API
accepts arbitrary complex domains.

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

Mixed-sign weight-two words may retain exact
`Li₂((1+i)/2)` or `Li₂((1-i)/2)`. No unproved minimal constant basis is
asserted. Repeated-letter words use powers divided by factorials; classical
higher words retain exact `Li_n(±i)`.

Tests independently integrate the defining one-dimensional path integral for
all convergent weight-two cases to 35 digits, using the explicit logarithmic
base case for `H_00`. Exact principal-log/Catalan identities, higher simple
words, and unsupported-input rejection are also checked. The numerical
comparisons supplement the branch and integral derivation.

## Mixed words at weights three and four

Split the straight path at `i/2`. The path composition formula sums over
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
strict nested sums. The returned constants support `.rewrite(Sum)`
for their exact defining sums and `.evalf()` for numerical evaluation.

The live Wolfram test compares all 121 words of weights zero through four
against HPL 2.0 at 35 decimal digits. The exact expressions can use a
different constant basis from the upstream reductions.
