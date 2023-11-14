def kmers(seq, k=4):
    n_kmers = len(seq) - k + 1
    for i in range(n_kmers):
        yield seq[i:(i + k)]

def kscore(seq, k=4):
    return len(set(kmers(seq))) / len(seq)
