# Native expansion pipeline

```mermaid
flowchart TD
    A[Exact SymPy input] --> B[Validate assumptions and endpoints]
    B --> C[Plan Laurent pole budget]
    C --> D[Select canonical native kernel]
    D --> E[Apply contiguous shifts]
    E --> F[Combine complete prefactors]
    F --> G[Truncate and attach provenance]
```

Canonical kernels use strict nested sums, colored sums or rationalized HPL
systems. Numerical evaluation is a separate API and is never silently used as
an exact expansion backend. Unsupported cases fail explicitly.
