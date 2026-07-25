import os
import pytest

from dnastore.formats import (
    FastaRecord, FastqRecord, write_fasta, read_fasta, write_fastq, read_fastq,
    uniform_quality_string, to_seqrecord,
)


def test_fasta_round_trip(tmp_path):
    records = [
        FastaRecord(id="s1", sequence="ACGT" * 30, description="test record"),
        FastaRecord(id="s2", sequence="TTTTGGGGCCCCAAAA", description=""),
    ]
    path = str(tmp_path / "out.fasta")
    write_fasta(records, path)
    parsed = read_fasta(path)
    assert [(r.id, r.sequence) for r in parsed] == [(r.id, r.sequence) for r in records]


def test_fasta_wraps_long_sequences(tmp_path):
    seq = "A" * 200
    write_fasta([FastaRecord(id="s1", sequence=seq)], str(tmp_path / "out.fasta"))
    with open(tmp_path / "out.fasta") as f:
        lines = f.read().splitlines()
    # header + wrapped sequence lines, none exceeding FASTA_LINE_WIDTH
    seq_lines = lines[1:]
    assert all(len(line) <= 70 for line in seq_lines)
    assert "".join(seq_lines) == seq


def test_fastq_round_trip(tmp_path):
    quality = uniform_quality_string(16)
    records = [FastqRecord(id="r1", sequence="ACGTACGTACGTACGT", quality=quality, description="sim")]
    path = str(tmp_path / "out.fastq")
    write_fastq(records, path)
    parsed = read_fastq(path)
    assert parsed[0].sequence == records[0].sequence
    assert parsed[0].quality == quality
    assert len(parsed[0].quality) == len(parsed[0].sequence)


def test_uniform_quality_string_length():
    q = uniform_quality_string(50, phred_score=30)
    assert len(q) == 50
    assert all(c == q[0] for c in q)


def test_biopython_interop_raises_clear_error_without_biopython():
    record = FastaRecord(id="s1", sequence="ACGT")
    with pytest.raises(ImportError, match="biopython"):
        to_seqrecord(record)
