# Exact Taylor series at an interior point

`hpl_taylor(word, center, variable, through=N)` returns exact coefficients
through degree $N$, plus an explicit SymPy `Order` term. The base point must
be an exact rational $0<c<1$. Coefficients retain HPL values at $c$ as exact
constants, rather than approximating them numerically.

Write $H_w(c+h)=\sum_n C_{w,n}h^n$ and
$f_a(c+h)=\sum_j K_{a,j}h^j$. The defining differential equation gives

$$
C_{w,0}=H_w(c),\qquad
C_{(a,v),n}=\frac1n\sum_{j=0}^{n-1}K_{a,j}C_{v,n-1-j}.
$$

The empty word has $C_{\varnothing,0}=1$ and zero higher coefficients. The
three kernel series are geometric:

$$
K_{0,j}=\frac{(-1)^j}{c^{j+1}},\qquad
K_{1,j}=\frac1{(1-c)^{j+1}},\qquad
K_{-1,j}=\frac{(-1)^j}{(1+c)^{j+1}}.
$$

The recurrence is finite and exact at every requested order. All three kernel
series, and hence the selected analytic HPL branch, converge in at least the
disk $|h|<\min(c,1-c)$. This is not an endpoint asymptotic expansion or a
prescription to cross a branch cut.

```python title="Expand around one third without numerical differentiation" linenums="1" hl_lines="5 6"
import sympy as s
from hypra import HPL, hpl_taylor

z, h = s.symbols("z h")
local = hpl_taylor((0, 1), s.Rational(1, 3), z, through=3)
composed = HPL((0, 1), s.Rational(1, 3) + h**2).series(h, 0, 5)
```

SymPy's `.series()` uses the same recurrence for arguments tending to a
supported interior base point with a known positive rational valuation.
The explicit API supports weight at most 24 and degree at most 64; unsupported
centers and resource requests raise errors rather than returning partial data.
