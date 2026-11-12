# Stuffle algebra of endpoint constants

`stuffle_expand_constants(expression)` linearizes exact polynomial expressions
in `MZV` and `EulerSum`. This is the **strict nested-sum product**, not the HPL
shuffle product and not a divergent-series regularization.

At a common finite cutoff, split the outer summation indices into the three
regions $n>m$, $m>n$ and $n=m$. The first two interleave the leading indices;
the third merges their weights and multiplies their colors. Recursing gives

$$
(a,A)*(b,B)=a\,(A*(b,B))+b\,((a,A)*B)+(a\circ b)\,(A*B),
$$

where a letter is a pair $a=(s,c)$ and
$(s,c)\circ(t,d)=(s+t,cd)$. Every term has the same total weight. Because
this identity is finite-cutoff exact, taking the convergent outer-cutoff limits
preserves it, including the allowed alternating sums with leading index one.
No undefined $\zeta(1)$ is introduced.

```python title="Linearize an exact zeta product" linenums="1" hl_lines="4"
from hypra import MZV, stuffle_expand_constants

zeta2 = MZV((2,))
result = stuffle_expand_constants(zeta2**2)
# 2*MZV((2, 2)) + MZV((4,))
```

A term with all colors $+1$ is returned as `MZV`; otherwise it remains an exact
`EulerSum`. Reducing those constants to ordinary zeta values is a separate,
optional operation. Arbitrary functions or inverse powers of endpoint constants
are rejected rather than treated as polynomial products. The implementation
limits total weight to twelve, powers to eight, and intermediate terms to 8192.
This is a linearization algorithm, not a claim of a minimal transcendental basis.
