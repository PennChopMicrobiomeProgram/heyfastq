from .read import Read


def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF
    args = [iter(iterable)] * n
    return zip(*args)


def parse_fastq(f):
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        yield Read(desc, seq, qual)


def parse_fastq_paired(fs):
    fastq_iters = [parse_fastq(f) for f in fs]
    for paired_rec in zip(*fastq_iters):
        yield paired_rec


def write_fastq(f, reads):
    for read in reads:
        f.write("@{0}\n{1}\n+\n{2}\n".format(read.desc, read.seq, read.qual))


def write_fastq_paired(fs, paired_recs):
    fs = list(fs)
    for paired_rec in paired_recs:
        for f, read in zip(fs, paired_rec):
            f.write("@{0}\n{1}\n+\n{2}\n".format(read.desc, read.seq, read.qual))


def parse_seq_ids(f):
    for line in f:
        line = line.strip()
        if line.startswith("#") or (line == ""):
            continue
        seq_id = line.split()[0]
        yield seq_id
