"""Demonstrate exporting a DNAStorage archive to FASTA (synthesis submission
format) and FASTQ (sequencer readback format) -- the two formats that let
this project's simulated pipeline talk to the rest of the bioinformatics
ecosystem (BLAST, alignment tools, real synthesis vendor order forms, etc).

Run: python examples/bio_interop.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dnastore import DNAStorage


def main():
    dna = DNAStorage(codec="rotating")
    original = os.urandom(3000)
    dna.store(original, name="dataset.bin")

    # FASTA: what you'd actually submit to a real synthesis vendor's order
    # form. Each record's description carries the addressing metadata
    # (primer pair, stripe/block index) that a real oligo pool would
    # encode via barcodes instead.
    n_fasta = dna.export_object_fasta("dataset.bin", "/tmp/dataset_synthesis_order.fasta")
    print(f"Exported {n_fasta} strands to FASTA (synthesis submission format)")

    # The whole physical pool, across every object -- useful for auditing
    # or feeding into general-purpose sequence analysis tools.
    n_pool = dna.export_pool_fasta("/tmp/entire_pool.fasta")
    print(f"Exported {n_pool} total strands across the whole archive")

    # FASTQ: what a real sequencer (or basecalling pipeline) would hand
    # back to you -- raw reads, post error-injection, with quality scores.
    # This is the format a downstream bioinformatics tool (aligners, QC
    # tools, custom analysis scripts) actually expects to consume.
    n_fastq = dna.retrieve_as_fastq("dataset.bin", "/tmp/dataset_reads.fastq")
    print(f"Exported {n_fastq} raw reads to FASTQ (sequencer output format)")

    print("\nFirst FASTA record:")
    with open("/tmp/dataset_synthesis_order.fasta") as f:
        print("".join(f.readlines()[:3]))

    print("With biopython installed (pip install dnastore[bio]), you can also do:")
    print("  from dnastore.formats import to_seqrecord, read_fasta")
    print("  records = [to_seqrecord(r) for r in read_fasta('/tmp/dataset_synthesis_order.fasta')]")
    print("  # now records are Bio.SeqRecord objects, usable with any Biopython-based tool")


if __name__ == "__main__":
    main()
