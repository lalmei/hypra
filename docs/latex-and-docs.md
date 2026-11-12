# Epsilon Expansions

Hypra expands the regulator $\varepsilon$ while keeping the argument and
coefficients symbolic. The generalized hypergeometric function is

$$
{}_pF_q(\mathbf a;\mathbf b;z)
=\sum_{n=0}^{\infty}\frac{\prod_i(a_i)_n}{\prod_j(b_j)_n}\frac{z^n}{n!}.
$$

For affine parameters $a_i=A_i+\alpha_i\varepsilon$, the implementation expands
Pochhammer coefficients, organizes harmonic sums, and converts their generating
functions to harmonic polylogarithms. Integer shifts use contiguous differential
operators; half-integer systems use rationalizing substitutions and differential
recurrences. It does not fit coefficients numerically.

## Example

$$
{}_2F_1(1,-\varepsilon;1-2\varepsilon;z)
=1+\varepsilon\log(1-z)
+\varepsilon^2\left[-2\operatorname{Li}_2(z)
-\frac12\log^2(1-z)\right]+O(\varepsilon^3).
$$

The convention is $H_1(z)=-\log(1-z)$ and $H_{0,1}(z)=\operatorname{Li}_2(z)$.

## Poles require deeper internal orders

If an external prefactor starts at $\varepsilon^{-r}$, a requested result
through $\varepsilon^N$ needs the regular factor through at least
$\varepsilon^{N+r}$. Hypra plans the complete expression before truncating.
It rejects unresolved orders, unsupported branches and unjustified endpoint
limits instead of returning a partial expansion as complete.

See [the detailed derivation](mathematics.md), [integer reductions](integer-port.md)
and [half-integer systems](half-system.md) for proofs, assumptions and worked
examples. [The PDF](mathematics.pdf) and LaTeX source are also in the repository.

## Local rendering

Equations use locally bundled KaTeX with LaTeX syntax. No browser CDN is needed.
