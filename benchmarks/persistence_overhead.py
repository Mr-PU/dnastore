"""
Store/retrieve wall-clock timing WITH on-disk persistence enabled
(state_path/pool_path set).

Unlike benchmarks/codec_comparison.py and benchmarks/efficiency_report.py,
which construct DNAStorage() with no path (in-memory only), this benchmark
measures the one thing those never touch: PhysicalPool's on-disk save
behavior on every store().

Mirrors the same (codec, size) cells published in the README's
performance table so results are comparable in shape -- though not
in scenario, since that table never enabled persistence at all.

Run: python benchmarks/persistence_overhead.py [label]
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dnastore import DNAStorage

CELLS = [
    ("fountain", "50 KB", 50 * 1024),
    ("fountain", "200 KB", 200 * 1024),
    ("fountain", "500 KB", 500 * 1024),
    ("rotating", "50 KB", 50 * 1024),
    ("rotating", "200 KB", 200 * 1024),
    ("naive", "500 KB", 500 * 1024),
]


def run_cell(codec_name: str, size_bytes: int) -> tuple[float, float, int]:
    tmpdir = tempfile.mkdtemp(prefix="dnastore_bench_")
    try:
        dna = DNAStorage(
            codec=codec_name,
            state_path=os.path.join(tmpdir, "state.json"),
            pool_path=os.path.join(tmpdir, "pool.json"),
        )
        data = os.urandom(size_bytes)

        t0 = time.perf_counter()
        info = dna.store(data, name="bench.bin")
        store_time = time.perf_counter() - t0

        t0 = time.perf_counter()
        out = dna.retrieve("bench.bin")
        retrieve_time = time.perf_counter() - t0

        assert out == data, f"round-trip mismatch: {codec_name} @ {size_bytes} bytes"
        return store_time, retrieve_time, info.strand_count
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    label = sys.argv[1] if len(sys.argv) > 1 else ""
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=os.path.dirname(__file__), text=True
        ).strip()
    except Exception:
        commit = "unknown"

    print(f"Persistence overhead benchmark{' -- ' + label if label else ''}")
    print(f"commit: {commit}   time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("state_path/pool_path ARE set for every run below (disk persistence ON) --")
    print("this is the scenario the README's own performance table never measures.")
    print()
    print(f"{'codec':<10} {'size':<8} {'strands':<9} {'store (s)':<12} {'retrieve (s)':<12}")
    print("-" * 55)
    sys.stdout.flush()

    for codec_name, size_label, size_bytes in CELLS:
        store_t, retrieve_t, n_strands = run_cell(codec_name, size_bytes)
        print(f"{codec_name:<10} {size_label:<8} {n_strands:<9} {store_t:<12.3f} {retrieve_t:<12.3f}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
