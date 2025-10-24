from typing import Dict, Optional


def kmers(seq, k=4):
    n_kmers = len(seq) - k + 1
    for i in range(n_kmers):
        yield seq[i : (i + k)]


_BASE_TO_BITS: Dict[str, int] = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "a": 0,
    "c": 1,
    "g": 2,
    "t": 3,
}


_BITSET_CACHE: Dict[int, bytearray] = {}
_MAX_BITSET_SIZE = 1 << 22  # 4 Mi entries (4 MiB) per cached bitset


def _kscore_fallback(seq: str, k: int) -> float:
    seq_len = len(seq)
    if seq_len == 0:
        return 0.0
    return len(set(kmers(seq, k=k))) / seq_len


def _bitset_for_k(k: int) -> Optional[bytearray]:
    size = 1 << (2 * k)
    if size > _MAX_BITSET_SIZE:
        return None

    bitset = _BITSET_CACHE.get(k)
    if bitset is None or len(bitset) != size:
        bitset = bytearray(size)
        _BITSET_CACHE[k] = bitset
    return bitset


def kscore(seq, k=4):
    seq_len = len(seq)
    if seq_len == 0:
        return 0.0

    if k <= 0:
        return _kscore_fallback(seq, k)

    if seq_len < k:
        return 0.0

    bitset = _bitset_for_k(k)
    if bitset is None:
        return _kscore_fallback(seq, k)

    mask = (1 << (2 * k)) - 1
    touched = []
    unique = 0
    rolling_value = 0
    window_len = 0
    use_fallback = False

    for base in seq:
        bits = _BASE_TO_BITS.get(base)
        if bits is None:
            use_fallback = True
            break

        rolling_value = ((rolling_value << 2) | bits) & mask
        window_len += 1

        if window_len >= k:
            idx = rolling_value
            if not bitset[idx]:
                bitset[idx] = 1
                touched.append(idx)
                unique += 1

    for idx in touched:
        bitset[idx] = 0

    if use_fallback:
        return _kscore_fallback(seq, k)

    return unique / seq_len
