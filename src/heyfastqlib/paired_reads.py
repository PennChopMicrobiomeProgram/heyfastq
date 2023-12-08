import operator

from .seqs import kscore

def trim_fixed(read, length=100):
    name, seq, qual = read
    return name, seq[:length], qual[:length]

def kscore_ok(read, k=4, min_kscore=0.55):
    name, seq, qual = read
    return kscore(seq, k=k) >= min_kscore

def read_length(read):
    return len(read[2])

def compare_length(read, threshold=100, cmp=operator.ge):
    return cmp(read_length(read), threshold)

def map_paired(paired_reads, f, **kwargs):
    for paired_read in paired_reads:
        yield tuple(f(read, **kwargs) for read in paired_read)

def filter_paired(paired_reads, f, requirement=all, **kwargs):
    for paired_read in paired_reads:
        passes_filter = (f(read, **kwargs) for read in paired_read)
        if requirement(passes_filter):
            yield paired_read
