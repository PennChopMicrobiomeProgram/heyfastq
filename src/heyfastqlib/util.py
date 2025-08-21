import collections
import itertools
import random


def subsample(xs, n, seed=None):
    random.seed(seed)
    # https://en.wikipedia.org/wiki/Reservoir_sampling
    reservoir = [None for _ in range(n)]
    for i, x in enumerate(xs):
        if i < n:
            reservoir[i] = x
        else:
            idx = random.randint(0, i)
            if idx < n:
                reservoir[idx] = x
    return reservoir


def sliding_sum(xs, k=4):
    # From moving_average recipe in Python docs
    xs = iter(xs)
    d = collections.deque(itertools.islice(xs, k - 1))
    d.appendleft(0)
    s = sum(d)
    for elem in xs:
        s += elem - d.popleft()
        d.append(elem)
        yield s
