"""Quickstart: basic store/retrieve/update/delete against a simulated
DNA storage backend. Run: python examples/quickstart.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dnastore import DNAStorage


def main():
    # `fountain` is the default codec: self-redundant, tolerant of
    # strand dropout with no separate ECC layer needed.
    dna = DNAStorage(
        codec="fountain",
        state_path="/tmp/dnastore_example_state.json",
        pool_path="/tmp/dnastore_example_pool.json",
    )

    # Store some in-memory data (store() also accepts a real file path,
    # e.g. dna.store("report.pdf"))
    original = b"Hello from DNA storage! " * 100
    info = dna.store(original, name="hello.txt")
    print(f"Stored '{info.name}' as {info.strand_count} DNA strands "
          f"({info.size_bytes} bytes, codec={info.codec})")

    # Retrieve it back -- this triggers simulated synthesis/sequencing
    # error injection under the hood, and the fountain decoder recovers
    # from whatever strand loss occurred.
    recovered = dna.retrieve("hello.txt")
    print("Round-trip successful:", recovered == original)

    # Update creates a new version; the old strands are marked obsolete
    # (not physically erased -- DNA can't be rewritten in place).
    dna.update(b"Updated content!", name="hello.txt")
    print("Versions:", dna.list_versions("hello.txt"))
    print("Latest content:", dna.retrieve("hello.txt"))

    # Delete tombstones the object; retrieve() will raise afterward.
    dna.delete("hello.txt")
    print("Objects after delete:", dna.list_objects())

    print("\nPool stats:", dna.stats())


if __name__ == "__main__":
    main()
