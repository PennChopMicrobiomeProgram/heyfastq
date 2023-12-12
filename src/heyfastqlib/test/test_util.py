import random

from heyfastqlib import util


def test_subsample():
    xs = ["a", "b", "c", "d", "e", "f"]
    # sampling less than or equal to n is original list
    assert util.subsample(xs[:4], n=4) == xs[:4]
    random.seed(0)
    # if any of the first n elements are in the resulting list,
    # they should be in their original position
    assert util.subsample(xs, n=3) == ["f", "b", "c"]


def test_sliding_sum():
    xs = [1, 2, 3, 4, 5]
    assert list(util.sliding_sum(xs, 3)) == [6, 9, 12]
