import collections
import itertools
import random
from typing import Generator, Optional


def subsample(xs: list[int], n: int, seed: Optional[int] = None) -> list[int]:
    random.seed(seed)
    # https://en.wikipedia.org/wiki/Reservoir_sampling
    reservoir: list[Optional[int]] = [None for _ in range(n)]
    for i, x in enumerate(xs):
        if i < n:
            reservoir[i] = x
        else:
            idx = random.randint(0, i)
            if idx < n:
                reservoir[idx] = x

    full_reservoir = [x for x in reservoir if x is not None]
    assert (
        len(full_reservoir) == n
    ), f"Expected reservoir of size {n}, got {len(full_reservoir)}"
    return full_reservoir


def sliding_sum(xs: list[int], k: int = 4) -> Generator[int, None, None]:
    # From moving_average recipe in Python docs
    it = iter(xs)
    d = collections.deque(itertools.islice(it, k - 1))
    d.appendleft(0)
    s = sum(d)
    for elem in it:
        s += elem - d.popleft()
        d.append(elem)
        yield s
