from dataclasses import dataclass
import operator
from typing import Callable, Generator, Tuple, TypeVar
from .seqs import kscore
from .util import sliding_sum


@dataclass(slots=True)
class Read:
    """A single FASTQ read representation"""

    desc: str
    seq: str
    qual: str


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


def seq_id(read: Read) -> str:
    return read.desc.split()[0]


def length(read: Read) -> int:
    return len(read.seq)


def qvals(read: Read, offset: int = 33) -> list[int]:
    return [ord(x) - offset for x in read.qual]


def trim_read(read: Read, start_idx: int, end_idx: int) -> Read:
    seq = read.seq[start_idx:end_idx]
    qual = read.qual[start_idx:end_idx]
    return Read(read.desc, seq, qual)


def trim(r: R, start_idx: int = 0, end_idx: int = 100) -> R:
    if isinstance(r, Read):
        return trim_read(r, start_idx, end_idx)
    elif isinstance(r, tuple):
        return (
            trim_read(r[0], start_idx, end_idx),
            trim_read(r[1], start_idx, end_idx),
        )


def kscore_ok(r: R, k: int = 4, min_kscore: float = 0.55) -> bool:
    if isinstance(r, Read):
        return kscore(r.seq, k=k) >= min_kscore
    elif isinstance(r, tuple):
        return all(kscore_ok(read, k=k, min_kscore=min_kscore) for read in r)


def length_ok(r: R, threshold: int = 100, cmp: Callable = operator.ge) -> bool:
    if isinstance(r, Read):
        return cmp(length(r), threshold)
    elif isinstance(r, tuple):
        return all(length_ok(read, threshold=threshold, cmp=cmp) for read in r)


def seq_id_ok(r: R, seq_ids: set[str] = set(), keep: bool = False) -> bool:
    if isinstance(r, Read):
        ids = set([seq_id(r)])
    elif isinstance(r, tuple):
        ids = set([seq_id(read) for read in r])

    if keep:
        return all(id in seq_ids for id in ids)
    else:
        return all(id not in seq_ids for id in ids)


def trim_read_moving_average(read: Read, k: int = 4, threshold: int = 15) -> Read:
    window_sum_threshold = threshold * k
    qs = qvals(read)
    for window_idx, window_sum in enumerate(sliding_sum(qs, k=k)):
        if window_sum < window_sum_threshold:
            end_idx = window_idx
            # Extend to include last qval in window meeting threshold
            for extended_idx in range(window_idx, window_idx + k):
                if qs[extended_idx] >= threshold:
                    end_idx = extended_idx + 1
            return trim(read, end_idx=end_idx)
    return read


def trim_moving_average(r: R, k: int = 4, threshold: int = 15) -> R:
    if isinstance(r, Read):
        return trim_read_moving_average(r, k=k, threshold=threshold)
    elif isinstance(r, tuple):
        return (
            trim_read_moving_average(r[0], k=k, threshold=threshold),
            trim_read_moving_average(r[1], k=k, threshold=threshold),
        )


def trim_read_ends(
    read: Read, threshold_start: int = 3, threshold_end: int = 3
) -> Read:
    qs = qvals(read)
    trim_start = 0
    for i, q in enumerate(qs):
        if q < threshold_start:
            trim_start = i + 1
        else:
            break
    trim_end = 0
    for i, q in enumerate(reversed(qs)):
        if q < threshold_end:
            trim_end = i + 1
        else:
            break
    if (trim_start == 0) and (trim_end == 0):
        return read
    else:
        start_idx = trim_start
        end_idx = length(read) - trim_end
        return trim(read, start_idx=start_idx, end_idx=end_idx)


def trim_ends(r: R, threshold_start: int = 3, threshold_end: int = 3) -> R:
    if isinstance(r, Read):
        return trim_read_ends(
            r, threshold_start=threshold_start, threshold_end=threshold_end
        )
    elif isinstance(r, tuple):
        return (
            trim_read_ends(
                r[0], threshold_start=threshold_start, threshold_end=threshold_end
            ),
            trim_read_ends(
                r[1], threshold_start=threshold_start, threshold_end=threshold_end
            ),
        )
