"""Standalone smoke test -- exercises the same scenarios as tests/ without
requiring pytest, so it can run in offline/no-network environments."""
import os
import random
import sys
import traceback

sys.path.insert(0, os.path.dirname(__file__))

from dnastore import DNAStorage
from dnastore.codec import NaiveCodec, RotatingCodec, FountainCodec, get_codec, available_codecs
from dnastore.ecc.reed_solomon import ReedSolomonErasureCoder
from dnastore.formats import (
    FastaRecord, FastqRecord, write_fasta, read_fasta, write_fastq, read_fastq,
    uniform_quality_string, to_seqrecord,
)
from dnastore.plugins import discover_plugin_codecs

passed = 0
failed = 0


def check(name, condition):
    global passed, failed
    if condition:
        print(f"PASS  {name}")
        passed += 1
    else:
        print(f"FAIL  {name}")
        failed += 1


def check_raises(name, exc_type, fn, *args, **kwargs):
    global passed, failed
    try:
        fn(*args, **kwargs)
        print(f"FAIL  {name} (no exception raised)")
        failed += 1
    except exc_type:
        print(f"PASS  {name}")
        passed += 1
    except Exception as e:
        print(f"FAIL  {name} (wrong exception: {e!r})")
        failed += 1


def main() -> int:
    global passed, failed
    passed = 0
    failed = 0
    # --- Codecs ---
    for codec in [NaiveCodec(), RotatingCodec(), FountainCodec()]:
        data = os.urandom(500)
        strands = codec.encode(data)
        decoded = codec.decode(strands, len(data))
        check(f"codec round-trip [{codec.name}]", decoded == data)

        empty_strands = codec.encode(b"")
        empty_decoded = codec.decode(empty_strands, 0)
        check(f"codec empty round-trip [{codec.name}]", empty_decoded == b"")

    rot = RotatingCodec()
    strands = rot.encode(bytes([0] * 1000))
    check("rotating has no homopolymers", all(rot.max_homopolymer_run(s) == 1 for s in strands))

    naive = NaiveCodec()
    strands = naive.encode(bytes([0] * 100))
    check("naive can produce homopolymers", any(naive.max_homopolymer_run(s) > 1 for s in strands))

    fountain = FountainCodec(redundancy=2.5)  # matches DNAStorage's default margin
    data = os.urandom(2000)
    strands = fountain.encode(data)
    surviving = [s for s in strands if random.random() > 0.1]
    decoded = fountain.decode(surviving, len(data))
    check("fountain tolerates 10% dropout", decoded == data)

    check_raises("get_codec unknown raises", ValueError, get_codec, "not-a-real-codec")

    # --- Reed-Solomon ECC ---
    data = os.urandom(1000)
    rs = ReedSolomonErasureCoder(data_shards=10, parity_shards=4)
    blocks = rs.encode(data)
    check("RS encode produces k+m blocks", len(blocks) == 14)
    check("RS decode no loss", rs.decode(blocks, len(data)) == data)

    random.shuffle(blocks)
    surviving = blocks[4:]
    check("RS decode tolerates max erasures", rs.decode(surviving, len(data)) == data)

    random.shuffle(blocks)
    check_raises("RS decode fails below threshold", ValueError, rs.decode, blocks[5:], len(data))

    # --- Full API ---
    for codec_name in ["naive", "rotating", "fountain"]:
        dna = DNAStorage(codec=codec_name)
        data = os.urandom(2000)
        dna.store(data, name="file.bin")
        out = dna.retrieve("file.bin")
        check(f"api store/retrieve [{codec_name}]", out == data)

        dna.update(b"v2", name="file.bin")
        versions = dna.list_versions("file.bin")
        check(f"api versioning [{codec_name}]", len(versions) == 2 and dna.retrieve("file.bin") == b"v2")

        dna.delete("file.bin")
        check(f"api delete tombstones [{codec_name}]", "file.bin" not in dna.list_objects())
        check_raises(f"api retrieve after delete raises [{codec_name}]", FileNotFoundError, dna.retrieve, "file.bin")

        stats_before = dna.stats()
        dna.store(b"more data", name="another.bin")
        dna.purge("another.bin")
        check(f"api purge removes strands [{codec_name}]", "another.bin" not in dna.list_objects(include_tombstoned=True))

    check_raises("retrieve missing object raises", FileNotFoundError, DNAStorage().retrieve, "nope.txt")
    check_raises("store bytes without name raises", ValueError, DNAStorage().store, b"data")

    # --- FASTA/FASTQ formats ---
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "out.fasta")
        records = [FastaRecord(id="s1", sequence="ACGT" * 30, description="test"),
                   FastaRecord(id="s2", sequence="TTTTGGGGCCCCAAAA")]
        write_fasta(records, fasta_path)
        parsed = read_fasta(fasta_path)
        check("fasta round-trip", [(r.id, r.sequence) for r in parsed] == [(r.id, r.sequence) for r in records])

        fastq_path = os.path.join(tmpdir, "out.fastq")
        quality = uniform_quality_string(16)
        fq_records = [FastqRecord(id="r1", sequence="ACGTACGTACGTACGT", quality=quality)]
        write_fastq(fq_records, fastq_path)
        parsed_fq = read_fastq(fastq_path)
        check("fastq round-trip", parsed_fq[0].sequence == fq_records[0].sequence and parsed_fq[0].quality == quality)

    check_raises("biopython interop raises without biopython", ImportError, to_seqrecord, FastaRecord(id="s1", sequence="ACGT"))

    # --- DNAStorage FASTA/FASTQ export ---
    with tempfile.TemporaryDirectory() as tmpdir:
        dna = DNAStorage(codec="rotating")
        dna.store(os.urandom(1000), name="sample.bin")

        n = dna.export_object_fasta("sample.bin", os.path.join(tmpdir, "obj.fasta"))
        check("export_object_fasta writes records", n > 0)

        n_pool = dna.export_pool_fasta(os.path.join(tmpdir, "pool.fasta"))
        check("export_pool_fasta writes records", n_pool >= n)

        n_fq = dna.retrieve_as_fastq("sample.bin", os.path.join(tmpdir, "reads.fastq"))
        check("retrieve_as_fastq writes reads", n_fq > 0)

        fastq_records = read_fastq(os.path.join(tmpdir, "reads.fastq"))
        check("fastq reads have matching seq/quality lengths",
              all(len(r.sequence) == len(r.quality) for r in fastq_records))

    # --- Plugin discovery ---
    import dnastore as dnastore_pkg

    codecs = dnastore_pkg.available_codecs()
    check("built-in codecs always present", {"naive", "rotating", "fountain"}.issubset(set(codecs)))
    check("discover_plugin_codecs returns a dict", isinstance(discover_plugin_codecs(), dict))
    check_raises("unknown codec raises with helpful message", ValueError, get_codec, "not-a-real-codec")

    if "xor_demo" in codecs:
        print("INFO  demo plugin package (dnastore-xor-demo) is installed -- testing it")
        plugin_codec = get_codec("xor_demo")
        data = os.urandom(500)
        strands = plugin_codec.encode(data)
        decoded = plugin_codec.decode(strands, len(data))
        check("plugin codec round-trip (direct)", decoded == data)

        dna = DNAStorage(codec="xor_demo")
        dna.store(data, name="f.bin")
        check("plugin codec round-trip (full API)", dna.retrieve("f.bin") == data)
    else:
        print("INFO  demo plugin package not installed -- skipping plugin-specific checks")
        print("INFO  (see examples/plugin_example/ to install it and exercise this path)")

    print(f"\n--- {passed} passed, {failed} failed ---")
    return failed


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
