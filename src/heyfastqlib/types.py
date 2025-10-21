from .read import Read
from typing import Generator, Protocol, TextIO, Tuple, TypeVar


ReadPair = Tuple[Read, Read]
R = TypeVar("R", Read, ReadPair)
"""The primary type variable representing either a single Read or a paired Read tuple.
This is what flows through the pipelines."""
ReadPipe = Generator[R, None, None]
"""The primary pipeline type representing a generator of Reads or Read pairs."""


def count_bases(r: R) -> int:
    if isinstance(r, tuple):
        return len(r[0].seq)
    else:
        return len(r.seq)
