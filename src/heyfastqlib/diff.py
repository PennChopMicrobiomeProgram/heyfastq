import itertools
from dataclasses import dataclass

from .read import seq_id

@dataclass
class FastqDiff:
    seq_id: str
    seq1: str
    seq2: str
    offset: int

    def format(self):
        out1, middle, out2 = zip(*self._format_zip())
        yield self.seq_id
        yield "  " + "".join(out1)
        yield "  " + "".join(middle)
        yield "  " + "".join(out2)

    def _format_zip(self):
        for c1, c2 in self._sequence_zip():
            middle = "|" if c1 == c2 else " "
            yield (c1, middle, c2)
        
    def _sequence_zip(self):
        seq1 = self.seq1
        if self.offset < 0:
            seq1 = (" " * -self.offset) + seq1
        seq2 = self.seq2
        if self.offset > 0:
            seq2 = (" " * self.offset) + seq2
        return itertools.zip_longest(seq1, seq2, fillvalue=" ")

    @classmethod
    def from_reads(cls, r1, r2):
        read_id = seq_id(r1)
        seq1 = r1.seq
        seq2 = r2.seq
        offset = diff_offset(seq1, seq2)
        return cls(read_id, seq1, seq2, offset)


def diff_offset(seq1, seq2):
    # offset of seq2 in seq1
    # negative if seq2 is extended on the left
    # no gaps in alignment
    # no overhangs split between seq1 and seq2
    len1 = len(seq1)
    len2 = len(seq2)
    if len1 == len2:
        if seq1 == seq2:
            return 0
        return best_hamming(seq1, seq2)
    if len1 > len2:
        find1 = seq1.find(seq2)
        if find1 >= 0:
            return find1
        return best_hamming(seq1, seq2)
    if len2 > len1:
        find2 = seq2.find(seq1)
        if find2 >= 0:
            #    AGCTCC
            # GGCAGCTCC find2 = 3, offset = -3
            return -find2
        return -best_hamming(seq2, seq1)


def idx_minimum(xs):
    xs = list(xs)
    return xs.index(min(xs))


def best_offset(long_seq, short_seq):
    ds = list(hamming_dists(long_seq, short_seq))
    return idx_minimum(ds)


def offset_dists(long_seq, short_seq):
    overlap_len = len(short_seq)
    end_idx = len(long_seq) - overlap_len + 1
    for idx in range(end_idx):
        segment = long_seq[idx:overlap_len]
        yield sum(c1 != c2 for c1, c2 in zip(segment, short_seq))
