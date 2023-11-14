import random

def subsample(xs, n):
    # Reservoir sampling avoids reading input data into memory
    # https://en.wikipedia.org/wiki/Reservoir_sampling
    reservoir = [None for _ in range(n)]
    for i, x in enumerate(xs):
        if i < n:
            reservoir[i] = x
        else:
            idx = random.randint(0, i)
            if idx < n:
                reservoir[idx] = seq
    return reservoir
