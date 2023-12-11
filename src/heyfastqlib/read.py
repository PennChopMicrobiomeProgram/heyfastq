from dataclasses import dataclass
import operator

from .seqs import kscore
from .util import sliding_sum


@dataclass(slots=True)
class Read:
    desc: str
    seq: str
    qual: str


def length(read):
    return len(read.seq)

def qvals(read, offset=33):
    return [ord(x) - offset for x in read.qual]

def trim(read, length=100):
    read.seq = read.seq[:length]
    read.qual = read.qual[:length]
    return read

def kscore_ok(read, k=4, min_kscore=0.55):
    return kscore(read.seq, k=k) >= min_kscore

def length_ok(read, threshold=100, cmp=operator.ge):
    return cmp(length(read), threshold)

def trim_moving_average(read, k=4, threshold=15):
    sum_threshold = threshold * k
    qs = qvals(read)
    for i, s in enumerate(sliding_sum(qs, k=k)):
        if s < sum_threshold:
            # trimmomatic also trims values from the end where the
            # per-base quality is less than the threshold
            return trim(read, i)
    return read
