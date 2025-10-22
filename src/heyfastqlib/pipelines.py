from itertools import islice
from multiprocessing import Pool
from typing import Callable, Iterable, Iterator

from .read import count_bases, R, ReadPipe


CounterDict = dict[str, int]


def _merge_counters(dest: CounterDict, src: CounterDict) -> None:
    for key, value in src.items():
        dest[key] = dest.get(key, 0) + value


def _chunk_reads(rs: Iterable[R], chunk_size: int) -> Iterator[list[R]]:
    iterator = iter(rs)
    while True:
        chunk = list(islice(iterator, chunk_size))
        if not chunk:
            break
        yield chunk


def _filter_worker(
    args: tuple[list[R], Callable[[R], bool], dict],
) -> tuple[list[R], CounterDict]:
    chunk, f, kwargs = args
    chunk_counter: CounterDict = {
        "input_reads": 0,
        "input_bases": 0,
        "output_reads": 0,
        "output_bases": 0,
    }
    out: list[R] = []
    for r in chunk:
        chunk_counter["input_reads"] += 1
        bases = count_bases(r)
        chunk_counter["input_bases"] += bases
        if f(r, **kwargs):
            chunk_counter["output_reads"] += 1
            chunk_counter["output_bases"] += bases
            out.append(r)
    return out, chunk_counter


def _map_worker(
    args: tuple[list[R], Callable[[R], R], dict],
) -> tuple[list[R], CounterDict]:
    chunk, f, kwargs = args
    chunk_counter: CounterDict = {
        "input_reads": 0,
        "input_bases": 0,
        "output_reads": 0,
        "output_bases": 0,
    }
    out: list[R] = []
    for r in chunk:
        chunk_counter["input_reads"] += 1
        input_bases = count_bases(r)
        chunk_counter["input_bases"] += input_bases
        mapped = f(r, **kwargs)
        chunk_counter["output_reads"] += 1
        chunk_counter["output_bases"] += count_bases(mapped)
        out.append(mapped)
    return out, chunk_counter


def filter_reads(
    rs: ReadPipe[R],
    f: Callable[[R], bool],
    counter: CounterDict,
    *,
    threads: int = 1,
    chunk_size: int = 1000,
    **kwargs,
) -> ReadPipe[R]:
    """
    Primary filtering pipeline
    """
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if chunk_size < 1:
        raise ValueError("chunk_size must be at least 1")

    # Using a multiprocessing Pool with a single worker is substantially slower
    # than iterating locally because of process start-up and pickle overhead,
    # so we keep this fast path for the common single-threaded case.
    if threads == 1:
        for r in rs:
            counter["input_reads"] += 1
            bases = count_bases(r)
            counter["input_bases"] += bases
            if f(r, **kwargs):
                counter["output_reads"] += 1
                counter["output_bases"] += bases
                yield r
        return

    task_iter = ((chunk, f, kwargs) for chunk in _chunk_reads(rs, chunk_size))
    with Pool(processes=threads) as pool:
        for out_chunk, chunk_counter in pool.imap(
            _filter_worker, task_iter, chunksize=1
        ):
            _merge_counters(counter, chunk_counter)
            for r in out_chunk:
                yield r


def map_reads(
    rs: ReadPipe[R],
    f: Callable[[R], R],
    counter: CounterDict,
    *,
    threads: int = 1,
    chunk_size: int = 1000,
    **kwargs,
) -> ReadPipe[R]:
    """
    Primary mapping pipeline
    """
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if chunk_size < 1:
        raise ValueError("chunk_size must be at least 1")

    # See note in filter_reads about the single-threaded fast path.
    if threads == 1:
        for r in rs:
            counter["input_reads"] += 1
            input_bases = count_bases(r)
            counter["input_bases"] += input_bases
            mapped = f(r, **kwargs)
            counter["output_reads"] += 1
            counter["output_bases"] += count_bases(mapped)
            yield mapped
        return

    task_iter = ((chunk, f, kwargs) for chunk in _chunk_reads(rs, chunk_size))
    with Pool(processes=threads) as pool:
        for out_chunk, chunk_counter in pool.imap(_map_worker, task_iter, chunksize=1):
            _merge_counters(counter, chunk_counter)
            for r in out_chunk:
                yield r


