"""Named HPL 2 constants from weakly ordered signed harmonic sums."""

import sympy as s

from ..symbolic import EulerSum
from .alternating import reduce_euler_sum


def _weak_sum(word):
    # Partition the weakly ordered simplex into strict blocks of equal indices.
    # A merged block adds exponents and multiplies colors.
    first, *rest = word
    partitions = [((abs(first), 1 if first > 0 else -1),)]
    for index in rest:
        power, color = abs(index), 1 if index > 0 else -1
        partitions = [
            block
            for prefix in partitions
            for block in (
                (*prefix, (power, color)),
                (*prefix[:-1], (prefix[-1][0]+power, prefix[-1][1]*color)),
            )
        ]
    return reduce_euler_sum(s.Add(*(
        EulerSum(*zip(*partition, strict=True)) for partition in partitions
    )))


HPLs6 = _weak_sum((-5, -1))
# These signs follow the executable numerical definitions and their values.
# The package's usage strings print different final signs for these two.
HPLs7a = _weak_sum((-5, 1, 1))
HPLs7b = _weak_sum((5, -1, -1))
HPLs8a = _weak_sum((5, 3))
HPLs8b = _weak_sum((-7, -1))
HPLs8c = _weak_sum((-5, -1, -1, -1))
HPLs8d = _weak_sum((-5, -1, 1, 1))
