"""
Codec comparison benchmark.

Produces the kind of table/plot that's actually worth putting in a
README for a project like this: information density, GC-content
stability, homopolymer risk, and empirical error tolerance for each
codec, all measured against the same random payload and the same
error injector.

Run: python benchmarks/codec_comparison.py
"""
from __future__ import annotations

import os
import random
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dnastore.codec import NaiveCodec, RotatingCodec, FountainCodec


def bits_per_base(codec_name: str, data: bytes, strands: list[str]) -> float:
    total_bases = sum(len(s) for s in strands)
    if total_bases == 0:
        return 0.0
    return (len(data) * 8) / total_bases


def gc_stats(strands: list[str]) -> tuple[float, float]:
    from dnastore.codec.base import Codec
    values = [Codec.gc_content(s) for s in strands if s]
    if not values:
        return 0.0, 0.0
    return statistics.mean(values), statistics.pstdev(values) if len(values) > 1 else 0.0


def max_homopolymer(strands: list[str]) -> int:
    from dnastore.codec.base import Codec
    return max((Codec.max_homopolymer_run(s) for s in strands), default=0)


def dropout_tolerance(codec, data: bytes, drop_fractions: list[float], trials: int = 15) -> dict[float, float]:
    results = {}
    for frac in drop_fractions:
        successes = 0
        for _ in range(trials):
            strands = codec.encode(data)
            surviving = [s for s in strands if random.random() > frac]
            try:
                decoded = codec.decode(surviving, len(data))
                if decoded == data:
                    successes += 1
            except Exception:
                pass
        results[frac] = successes / trials
    return results


def main():
    random.seed(42)
    data = os.urandom(5000)

    codecs = {
        "naive": NaiveCodec(),
        "rotating": RotatingCodec(),
        "fountain": FountainCodec(redundancy=2.5),
    }

    print(f"{'codec':<10} {'bits/base':<10} {'n_strands':<10} {'gc_mean':<9} {'gc_std':<8} {'max_homopolymer':<16}")
    print("-" * 70)
    encoded = {}
    for name, codec in codecs.items():
        strands = codec.encode(data)
        encoded[name] = strands
        bpb = bits_per_base(name, data, strands)
        gc_mean, gc_std = gc_stats(strands)
        homopolymer = max_homopolymer(strands)
        print(f"{name:<10} {bpb:<10.2f} {len(strands):<10} {gc_mean:<9.2f} {gc_std:<8.3f} {homopolymer:<16}")

    print()
    print("Dropout tolerance (fraction of trials that decoded successfully):")
    print(f"{'codec':<10} " + " ".join(f"{int(f*100):>4}%" for f in [0.05, 0.1, 0.15, 0.2, 0.3]))
    for name, codec in codecs.items():
        # naive/rotating alone have no redundancy -- this benchmark measures
        # the *codec's own* tolerance, not the full ECC-wrapped pipeline
        # (see the DNAStorage API for the RS-erasure-wrapped naive/rotating
        # numbers, which are much better).
        results = dropout_tolerance(codec, data, [0.05, 0.1, 0.15, 0.2, 0.3])
        row = " ".join(f"{v*100:>4.0f}%" for v in results.values())
        print(f"{name:<10} {row}")

    print()
    print("Note: naive/rotating show ~0% tolerance here because they have no")
    print("built-in redundancy -- that's the whole reason the ECC layer exists.")
    print("Use DNAStorage(codec='naive'/'rotating') to see the full RS-wrapped")
    print("pipeline's actual tolerance, which is governed by parity_fraction.")


if __name__ == "__main__":
    main()
