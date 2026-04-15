import itertools
from dataclasses import dataclass

from .read import seq_id


@dataclass
class FastqDiff:
    seq_id: str
    seq1: str
    seq2: str
    # Coordinate on seq1 where seq2 begins
    offset: int

    @property
    def overlap_length(self):
        return min(len(self.seq1), len(self.seq2))

    # AGCTCGAT
    #   CTCG
    # offset = 2
    # left_trim = seq1[:2]
    # aligned1 = seq1[2:6]
    # right_trim = seq1[6:]

    @property
    def start_idx1(self):
        return max(self.offset, 0)

    @property
    def end_idx1(self):
        return self.start_idx1 + self.overlap_length

    @property
    def left_trim(self):
        return self.seq1[: self.start_idx1]

    @property
    def aligned1(self):
        return self.seq1[self.start_idx1 : self.end_idx1]

    @property
    def right_trim(self):
        return self.seq1[self.end_idx1 :]

    #   CTCG
    # AGCTCGAT
    # offset = -2
    # left_extend = seq2[:2]
    # overlap = seq2[2:6]
    # right_extend = seq2[6:]

    @property
    def start_idx2(self):
        return max(-self.offset, 0)

    @property
    def end_idx2(self):
        return self.start_idx2 + self.overlap_length

    @property
    def left_extend(self):
        return self.seq2[: self.start_idx2]

    @property
    def aligned2(self):
        return self.seq2[self.start_idx2 : self.end_idx2]

    @property
    def right_extend(self):
        return self.seq2[self.end_idx2 :]

    @property
    def mismatches(self):
        idx = self.start_idx1
        for c1, c2 in zip(self.aligned1, self.aligned2):
            if c1 != c2:
                yield (idx, c1, c2)
            idx += 1

    @property
    def mismatches_format(self):
        mm = [",".join(map(str, cs)) for cs in self.mismatches]
        return ";".join(mm)

    field_names = [
        "read_id",
        "left_trim",
        "left_extend",
        "right_trim",
        "right_extend",
        "mismatches",
    ]

    field_attrs = [
        "seq_id",
        "left_trim",
        "left_extend",
        "right_trim",
        "right_extend",
        "mismatches_format",
    ]

    def format_tsv(self):
        vals = [getattr(self, a) for a in self.field_attrs]
        return "\t".join(vals) + "\n"

    def format_alignment(self, indent=2):
        indent_str = " " * indent
        aligned1, bars, aligned2 = self._alignment_seqs()
        return (
            f"{self.seq_id}\n"
            f"{indent_str}{aligned1}\n"
            f"{indent_str}{bars}\n"
            f"{indent_str}{aligned2}\n"
        )

    def _alignment_seqs(self):
        a1, bars, a2 = zip(*self._alignment_zip())
        return ("".join(a1), "".join(bars), "".join(a2))

    def _alignment_zip(self):
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
        return best_offset(seq1, seq2)
    if len1 > len2:
        find1 = seq1.find(seq2)
        if find1 >= 0:
            return find1
        return best_offset(seq1, seq2)
    if len2 > len1:
        find2 = seq2.find(seq1)
        if find2 >= 0:
            #    AGCTCC
            # GGCAGCTCC find2 = 3, offset = -3
            return -find2
        return -best_offset(seq2, seq1)


def idx_minimum(xs):
    xs = list(xs)
    return xs.index(min(xs))


def best_offset(long_seq, short_seq):
    ds = list(offset_dists(long_seq, short_seq))
    return idx_minimum(ds)


def offset_dists(long_seq, short_seq):
    overlap_len = len(short_seq)
    end_idx = len(long_seq) - overlap_len + 1
    for idx in range(end_idx):
        segment = long_seq[idx:overlap_len]
        yield sum(c1 != c2 for c1, c2 in zip(segment, short_seq))
