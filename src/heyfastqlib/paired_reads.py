from .seqs import kscore

def trim_fixed(read, length=100):
    name, seq, qual = read
    return name, seq[:length], qual[:length]

def kscore_ok(read, min_kscore=0.55):
    name, seq, qual = read
    return kscore(seq) >= min_kscore

def map_paired(paired_reads, f, **kwargs):
    for paired_read in paired_reads:
        yield (f(read, **kwargs) for read in paired_read)

def filter_paired(paired_reads, f, **kwargs):
    for paired_read in paired_reads:
        passes_filter = (f(read, **kwargs) for read in paired_read)
        if all(passes_filter):
            yield paired_read
