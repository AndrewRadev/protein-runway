"""
General-purpose utilities that can be used for various purposes across the
project.
"""

import itertools

def batched(iterable, n, *, strict=False):
    """
    Batch elements from the iterable into tuples of length n.

    Taken from itertools in python 3.12:
    <https://docs.python.org/3/library/itertools.html#itertools.batched>

    This function can be replaced with itertools.batched in the event of a
    python update.
    """
    if n < 1:
        raise ValueError('n must be at least one')

    iterator = iter(iterable)

    while batch := tuple(itertools.islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch
