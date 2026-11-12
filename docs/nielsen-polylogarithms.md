# Nielsen generalized polylogarithms

`NielsenPolylog(n, p, z)` is an exact, lossless notation layer with positive
integer indices and the explicit convention

\[
S_{n,p}(z)=H_{\underbrace{0,\ldots,0}_{n},\underbrace{1,\ldots,1}_{p}}(z),
\qquad n,p\ge1.
\]

Its HPL weight is `n+p`, not `n`. In particular `S_{n,1}(z)=Li_{n+1}(z)`.
Construction retains the Nielsen symbol; it does not claim a general
classical-polylogarithm reduction or silently change branches.

```python
import sympy as s
from hypra.nielsen import (
    NielsenPolylog, nielsen_to_hpl, hpl_to_nielsen, reduce_nielsen,
)
from hypra import HPL
z = s.Symbol("z")
f = NielsenPolylog(2, 3, z)
f.as_hpl()             # HPL((0,0,1,1,1), z)
f.rewrite(HPL)         # same exact conversion
hpl_to_nielsen(f.as_hpl()) == f
reduce_nielsen(NielsenPolylog(2, 1, z))  # polylog(3, z)
s.diff(f, z)           # NielsenPolylog(1,3,z)/z
s.latex(f)             # S_{2,3}(z)
```

Both conversion functions operate on symbolic expressions. Reverse conversion
changes only expanded words with a nonempty zero block followed by a nonempty
one block. Other HPLs, including empty words, pure logs, and alternating-letter
words, remain untouched. Exact JSON AST serialization uses a `NielsenPolylog`
head with the two integer indices and the argument. Invalid or floating-point
indices are rejected before constructor caching. The native serializer supports
this head; the optional Wolfram bridge does not currently decode it.

## Integral normalization

On the real interval `0<z<1`, the convention is equivalently

\[
S_{n,p}(z)=\frac{(-1)^{n+p-1}}{(n-1)!p!}
\int_0^1\frac{dt}{t}\log^{n-1}(t)\log^p(1-zt).
\]

To see the normalization, the all-one HPL is
`(-log(1-u))^p/p!`. Integrating `n` successive zero kernels gives the integral
of that function against `log(z/u)^(n-1)/((n-1)!u)`. Substituting `u=zt`
gives the displayed signs and factorials. The real integrals converge for
positive indices; extension outside this interval uses the matching principal
HPL continuation, not a fresh implicit branch prescription.

The derivative is

\[
\frac{d}{dz}S_{n,p}(z)=
\begin{cases}
S_{n-1,p}(z)/z,&n>1,\\
(-\log(1-z))^p/(p!z),&n=1.
\end{cases}
\]

This is a generic-point expression; taking its removable limit at zero is a
separate operation. Only the `p=1` classical reduction is implemented here.
Notation conversion has a weight-10000 allocation guard; lower numerical or
algebraic HPL resource limits still apply when those engines are used.
Tests compare the derivative with the HPL kernel and independently integrate
the Nielsen integral for several `(n,p)`, rather than merely reversing the
conversion map. No upstream source or special-value table is copied.
