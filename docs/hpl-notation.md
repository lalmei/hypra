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
