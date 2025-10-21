from .types import count_bases, R, ReadPipe
from .util import subsample
from typing import Callable, Optional


def filter_reads(
    rs: ReadPipe[R], f: Callable[[R], bool], counter: dict[str, int], **kwargs
) -> ReadPipe[R]:
    """
    Primary filtering pipeline
    """
    for r in rs:
        counter["input_reads"] += 1
        counter["input_bases"] += count_bases(r)
        if f(r, **kwargs):
            counter["output_reads"] += 1
            counter["output_bases"] += count_bases(r)
            yield r


def map_reads(
    rs: ReadPipe[R], f: Callable[[R], R], counter: dict[str, int], **kwargs
) -> ReadPipe[R]:
    """
    Primary mapping pipeline
    """
    for r in rs:
        counter["input_reads"] += 1
        counter["input_bases"] += count_bases(r)
        yield f(r, **kwargs)

    counter["output_reads"] = counter["input_reads"]
    counter["output_bases"] = counter["input_bases"]


def subsample_reads(
    rs: ReadPipe[R], l: int, n: int, seed: Optional[int]
) -> ReadPipe[R]:
    """
    Subsample is a weird one because it requires working with a full list of reads' indexes in memory
    Doesn't really fit with the pipeline model
    """
    idxs = subsample(list(range(l)), n, seed)
    for i, r in enumerate(rs):
        if i in idxs:
            yield r
