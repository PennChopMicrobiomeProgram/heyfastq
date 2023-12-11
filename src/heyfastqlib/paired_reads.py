def map_paired(paired_reads, f, **kwargs):
    for paired_read in paired_reads:
        yield tuple(f(read, **kwargs) for read in paired_read)


def filter_paired(paired_reads, f, requirement=all, **kwargs):
    for paired_read in paired_reads:
        passes_filter = (f(read, **kwargs) for read in paired_read)
        if requirement(passes_filter):
            yield paired_read
