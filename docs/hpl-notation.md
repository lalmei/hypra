# Index and kernel notation

`expand_indices` and `compress_indices` convert between expanded letters and
signed compressed integer indices. A zero run preceding ±1 is compressed
into the corresponding signed run length; trailing zeros remain standalone.
For example `(0,1,0,0,-1,0)` compresses to `(2,-3,0)`.

`hpl_pm(word,z)` accepts discrete labels `"plus"`, `"minus"`, and zero, with
kernels f_plus=f_1+f_-1=2/(1-z²) and
f_minus=f_1-f_-1=2z/(1-z²). It returns an exact sum of ordinary HPLs, so no
special serialization or hidden evaluation convention is needed.
`pm_to_letters` and `letters_to_pm` return exact word-coefficient dictionaries
for this multilinear basis change and its inverse. The implementation expands
kernel tensor products, not copied upstream conversion rules. Weight bounds
limit the exponential output size; round-trip tests cover all ternary words
through weight four.

Use `HPLPM(("plus", 0, "minus"), z)` to keep the kernel labels in a
symbolic expression. Calling `.rewrite(HPL)` converts it to ordinary HPLs.

## HPL 2.0 minimal sets

`reduce_hpl_to_minimal_set(expression)` reduces ordinary `HPL` and
`HPLPM` atoms in their respective alphabets through weight six. For each
weight, the reducer forms shuffle equations from all pairs of smaller words.
It eliminates words with leading one or trailing zero first, followed by
the remaining words in lexicographic order. In the plus/minus alphabet,
minus takes the place of −1 and plus takes the place of 1. Solving the
equations gives the upstream basis choices without loading reduction tables.
The live comparison checks every polynomial in both bases against HPL 2.0.
See [Maître, section 6](https://arxiv.org/abs/hep-ph/0507152).
