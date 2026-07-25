"""
FASTA/FASTQ export and Biopython interoperability.

FASTA and FASTQ are the actual lingua franca of DNA storage and
bioinformatics generally: every real synthesis vendor accepts FASTA
submissions, and every sequencer (or basecalling pipeline) emits
FASTQ. Without this module, dnastore's simulated strand pool is
locked inside its own JSON format and can't talk to anything else in
the bioinformatics ecosystem -- this module is what makes it
interoperable.

Two independent layers:
  - Pure Python FASTA/FASTQ read/write (no dependencies -- always
    available).
  - Optional Biopython SeqRecord conversion, for code that already
    works in terms of Bio.SeqRecord/Bio.Seq (BLAST wrappers, alignment
    tools, etc.). Requires `pip install dnastore[bio]`.
"""
from __future__ import annotations

from dataclasses import dataclass

FASTA_LINE_WIDTH = 70  # conventional wrap width for FASTA sequence lines


@dataclass
class FastaRecord:
    id: str
    sequence: str
    description: str = ""


@dataclass
class FastqRecord:
    id: str
    sequence: str
    quality: str  # Phred+33 ASCII-encoded, same length as sequence
    description: str = ""


# ---------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------

def write_fasta(records: list[FastaRecord], path: str) -> None:
    with open(path, "w") as f:
        for rec in records:
            header = f">{rec.id} {rec.description}".rstrip()
            f.write(header + "\n")
            for i in range(0, len(rec.sequence), FASTA_LINE_WIDTH):
                f.write(rec.sequence[i:i + FASTA_LINE_WIDTH] + "\n")


def fasta_string(records: list[FastaRecord]) -> str:
    lines = []
    for rec in records:
        header = f">{rec.id} {rec.description}".rstrip()
        lines.append(header)
        for i in range(0, len(rec.sequence), FASTA_LINE_WIDTH):
            lines.append(rec.sequence[i:i + FASTA_LINE_WIDTH])
    return "\n".join(lines) + ("\n" if lines else "")


def read_fasta(path: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    current_id = None
    current_desc = ""
    current_seq: list[str] = []

    def flush():
        if current_id is not None:
            records.append(FastaRecord(id=current_id, sequence="".join(current_seq), description=current_desc))

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
                parts = header.split(None, 1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                current_seq = []
            else:
                current_seq.append(line.strip())
        flush()

    return records


# ---------------------------------------------------------------------
# FASTQ
# ---------------------------------------------------------------------

def write_fastq(records: list[FastqRecord], path: str) -> None:
    with open(path, "w") as f:
        for rec in records:
            header = f"@{rec.id} {rec.description}".rstrip()
            f.write(header + "\n")
            f.write(rec.sequence + "\n")
            f.write("+\n")
            f.write(rec.quality + "\n")


def read_fastq(path: str) -> list[FastqRecord]:
    records: list[FastqRecord] = []
    with open(path, "r") as f:
        while True:
            header_line = f.readline()
            if not header_line:
                break
            seq_line = f.readline().rstrip("\n")
            plus_line = f.readline()  # noqa: F841 -- the '+' separator, unused
            qual_line = f.readline().rstrip("\n")

            header = header_line.rstrip("\n")[1:]
            parts = header.split(None, 1)
            rec_id = parts[0]
            desc = parts[1] if len(parts) > 1 else ""
            records.append(FastqRecord(id=rec_id, sequence=seq_line, quality=qual_line, description=desc))
    return records


def uniform_quality_string(length: int, phred_score: int = 37) -> str:
    """Generate a flat-confidence Phred+33 quality string.

    This is a deliberate simplification, not a claim of realism: the
    simulator doesn't currently track *which* bases within a read were
    altered by substitution/indel injection, so it can't report
    per-base confidence that actually reflects where errors occurred.
    A real sequencer's quality scores vary a lot per base and
    correlate with actual error likelihood -- wiring that through from
    SequencingSimulator's error injection is a good next contribution
    (see CONTRIBUTING.md).
    """
    char = chr(33 + phred_score)
    return char * length


# ---------------------------------------------------------------------
# Biopython interop (optional)
# ---------------------------------------------------------------------

def to_seqrecord(record: FastaRecord):
    try:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError as e:
        raise ImportError(
            "Biopython interop requires the optional 'biopython' dependency. "
            "Install with: pip install dnastore[bio]"
        ) from e
    return SeqRecord(Seq(record.sequence), id=record.id, description=record.description)


def from_seqrecord(seqrecord) -> FastaRecord:
    return FastaRecord(id=seqrecord.id, sequence=str(seqrecord.seq), description=seqrecord.description)


def write_records_via_biopython(records: list[FastaRecord], path: str, fmt: str = "fasta") -> None:
    """Write via Biopython's SeqIO instead of this module's own writer --
    useful mainly for formats beyond plain FASTA/FASTQ that Biopython
    supports (e.g. GenBank) and that this module doesn't implement
    natively."""
    try:
        from Bio import SeqIO
    except ImportError as e:
        raise ImportError(
            "Biopython interop requires the optional 'biopython' dependency. "
            "Install with: pip install dnastore[bio]"
        ) from e
    seqrecords = [to_seqrecord(r) for r in records]
    SeqIO.write(seqrecords, path, fmt)
