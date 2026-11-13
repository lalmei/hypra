# Full near-one power/log expansions

`hpl_series_at_one(word, t, through=N)` expands `H_word(1-t)` for the explicit
real approach `t -> 0+`. It includes all powers through `t**N`, with exact
polynomials in `T=-log(t)` at each power. Unlike the endpoint-polynomial API,
it retains the subleading powers as well as the leading logarithms.

```python
import sympy as s
from hypra.hpl_series_one import hpl_series_at_one

t = s.Symbol("t", positive=True)
result = hpl_series_at_one((0, 1), t, through=2)
# pi**2/6 + t*(log(t)-1) + t**2*(log(t)/2-1/4)
# + O(t**3*log(t)**2)
```

The log-enhanced remainder is intentional: a plain `O(t**(N+1))` would not
bound the omitted logarithms. The returned conservative bound is
`O(t**(N+1)*log(t)**weight)`. Empty words are exactly one. The function is not
a generic limit routine or an automatic complex continuation.

## Exact recurrence

Write

\[
H_w(1-t)=\sum_{n\ge0}t^n C_{w,n}(T),\qquad T=-\log t.
\]

`C_{w,0}` is the previously defined exact shuffle endpoint polynomial. For
`n>=1`, differentiating `t^n C(T)` gives
`t^(n-1)*(n*C-C')`. If the first letter is `a` and its suffix is `v`, the
kernel ODE therefore gives `(n-D_T) C_{w,n}=R_{w,n}`, with

\[
R_{w,n}=\begin{cases}
-C_{v,n},&a=1,\\
-\sum_{m=0}^{n-1}C_{v,m},&a=0,\\
-\sum_{m=0}^{n-1}2^{-(n-m)}C_{v,m},&a=-1.
\end{cases}
\]

These are respectively the transformed kernels `-1/t`, `-1/(1-t)` and
`-1/(2-t)`. On finite polynomials the unique inverse is

\[
C_{w,n}=\sum_{k\ge0}\frac{D_T^k R_{w,n}}{n^{k+1}},
\]

where the sum terminates when the derivative vanishes. Bottom-up suffix
recursion initializes the empty word with constant one. No unknown
integration constant is discarded: every zero-power constant comes from the
exact endpoint polynomial.

The rational kernels have convergent power series near `t=0`; the only
singular kernel is the explicit simple pole `1/t`. Induction on word length
therefore gives convergent local power/log expansions with log degree at
most the word weight. Omitting powers beyond `N` produces the stated
log-enhanced bound for real positive `t`. This derivation does not exchange
a regulator limit with an unrelated improper integral.

## Deliberate resource and branch limits

- Expanded alphabet `0,1,-1`, weight at most eight.
- Inclusive power order from zero through 24.
- Coefficient-operation budget of 12000; oversized expressions fail closed.
- Exact EulerSum/MZV endpoint constants may remain symbolic.
- `t` is an ordinary commutative SymPy symbol, with positive approach specified
  by the API contract even if the symbol has no positivity assumption.
- No automatic replacement of `HPL.series`, no other approach direction,
  and no arbitrary argument composition are introduced by this API.

Tests compare elementary logarithms, an independent dilogarithm reflection
identity, all nine weight-two ternary differential recurrences, and mixed-word
high-precision numerical checks. Numerical checks are secondary evidence;
the coefficients are derived by the exact recurrence above.

## Combined finite limits with rational poles

`hpl_limit_at_one(expression, z)` now also handles polynomial combinations
of ternary HPLs with numeric rational-function coefficients. It computes the
maximum coefficient pole order $P$, expands every HPL through $t^P$, and only
then combines Laurent powers and logarithms. The neglected product remainder
is at worst $O(t\log^k t)$ after multiplication by a pole, hence vanishes.
Every negative power or nonconstant logarithmic coefficient must cancel
exactly before a finite result is returned.

```python title="Cancel a rational pole before taking the endpoint" linenums="1" hl_lines="5"
import sympy as s
from hypra import HPL, hpl_limit_at_one

z = s.Symbol("z", positive=True)
result = hpl_limit_at_one((HPL((-1,), z) - s.log(2)) / (1-z), z)
assert result == -s.Rational(1, 2)
```

The ternary limit path bounds individual HPL weight by eight, total monomial
weight by twelve, and rational pole order by sixteen. Unknown symbolic
coefficient denominators, inverse HPLs and unproved cancellations fail closed.
It does not assign regularized constants to divergent expressions.
