from dataclasses import dataclass
import operator

from .seqs import kscore
from .util import sliding_sum


@dataclass(slots=True)
class Read:
    desc: str
    seq: str
    qual: str


def seq_id(read):
    return read.desc.split()[0]


def length(read):
    return len(read.seq)


def qvals(read, offset=33):
    return [ord(x) - offset for x in read.qual]


def trim(read, start_idx=0, end_idx=100):
    seq = read.seq[start_idx:end_idx]
    qual = read.qual[start_idx:end_idx]
    return Read(read.desc, seq, qual)


def kscore_ok(read, k=4, min_kscore=0.55):
    return kscore(read.seq, k=k) >= min_kscore


def length_ok(read, threshold=100, cmp=operator.ge):
    return cmp(length(read), threshold)


def seq_id_ok(read, seq_ids, remove=False):
    id = seq_id(read)
    if remove:
        return id not in seq_ids
    else:
        return id in seq_ids


def trim_moving_average(read, k=4, threshold=15):
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


def trim_ends(read, threshold_start=3, threshold_end=3):
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
