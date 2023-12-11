from heyfastqlib.read import Read, length, length_ok
from heyfastqlib.paired_reads import *

reads = [
    (Read("a", "AAAG", "1234"), Read("a", "CTTT", "1234")),
    (Read("b", "ACGTACGT", "12345678"), Read("b", "ACG", "123")),
]

def test_map_paired():
    assert list(map_paired(reads, length)) == [(4, 4), (8, 3)]

def test_filter_paired():
    assert list(filter_paired(reads, length_ok, threshold = 4)) == \
        reads[:1]

    assert list(
        filter_paired(reads, length_ok, requirement=any, threshold = 4)
    ) == reads
