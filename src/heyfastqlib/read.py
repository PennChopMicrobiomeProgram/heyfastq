from dataclasses import dataclass
import operator

from .seqs import kscore


@dataclass(slots=True)
class Read:
    desc: str
    seq: str
    qual: str


def length(read):
    return len(read.seq)


def trim(read, length=100):
    read.seq = read.seq[:length]
    read.qual = read.qual[:length]
    return read


def kscore_ok(read, k=4, min_kscore=0.55):
    return kscore(read.seq, k=k) >= min_kscore


def length_ok(read, threshold=100, cmp=operator.ge):
    return cmp(length(read), threshold)
