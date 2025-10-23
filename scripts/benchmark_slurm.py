#!/usr/bin/env python3
"""Submit heyfastq benchmarking jobs to Slurm."""

from __future__ import annotations

import itertools
import re
import shlex
import subprocess
from pathlib import Path

INPUTS = [
    Path("/path/to/input_R1.fastq.gz"),
    Path("/path/to/input_R2.fastq.gz"),
]
CMD = "trim-length"
ARGS = "--length 100"
THREADS = [1, 4, 16, 32]
COMPRESSION = [1, 9]
MEMORY = "8G"

def _output_path(input_path: Path) -> Path:
    return input_path.with_name(f"OUT_{input_path.name}")


def _wrap_command(inputs: list[Path], outputs: list[Path], threads: int) -> str:
    command_parts = [
        "heyfastq",
        CMD,
        *shlex.split(ARGS),
        "--threads",
        str(threads),
        "--input",
        *(str(path) for path in inputs),
        "--output",
        *(str(path) for path in outputs),
    ]
    return " ".join(shlex.quote(part) for part in command_parts)


def _job_name(threads: int, compression: int) -> str:
    base = INPUTS[0].name.replace(".", "-")
    return f"{CMD}-t{threads}-c{compression}-{base}"


def submit_jobs() -> None:
    outputs = [_output_path(path) for path in INPUTS]
    job_ids = []

    for threads, compression in itertools.product(THREADS, COMPRESSION):
        wrap_command = _wrap_command(INPUTS, outputs, threads)
        job_name = _job_name(threads, compression)
        sbatch_command = [
            "sbatch",
            f"--job-name={job_name}",
            f"--cpus-per-task={threads}",
            f"--mem={MEMORY}",
            f"--export=ALL,HFQ_GZIP_COMPRESSION={compression}",
            "--wrap",
            wrap_command,
        ]
        print("Submitting:", " ".join(shlex.quote(arg) for arg in sbatch_command))
        result = subprocess.run(
            sbatch_command,
            check=True,
            capture_output=True,
            text=True,
        )
        stdout = result.stdout.strip()
        match = re.search(r"Submitted batch job (\d+)", stdout)
        if match:
            job_ids.append(match.group(1))
        else:
            print("Warning: unable to parse job id from sbatch output:", stdout)
    print("Submitted job IDs:", ", ".join(job_ids) if job_ids else "(none)")


def main() -> None:
    submit_jobs()


if __name__ == "__main__":
    main()
