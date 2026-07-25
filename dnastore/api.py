"""
DNAStorage: the developer-facing API.

    dna = DNAStorage()
    dna.store("report.pdf")
    data = dna.retrieve("report.pdf")
    dna.update("report.pdf")     # new version; old strands marked obsolete
    dna.delete("report.pdf")     # tombstone; DNA is not physically erased
    dna.purge("report.pdf")     # explicit physical erase of all versions

Design notes (see README for the full rationale):
- Non-fountain codecs (naive/rotating) get their redundancy from a
  Reed-Solomon erasure code applied across strands *before* per-strand
  DNA encoding: the original bytes are split into k data shards + m
  parity shards, and each shard is encoded into exactly one DNA
  strand.
- The fountain codec is rateless/self-redundant and skips the RS
  layer entirely.
- Every strand gets a CRC16 appended (naive-encoded, 8 extra bases).
  A strand whose checksum doesn't match on readback is treated as
  *erased*, not corrected in place -- correction is delegated
  entirely to the cross-strand ECC layer (or the fountain decoder's
  own redundancy), which is a much better model of how real DNA
  storage systems actually recover from errors.
- Strand identity (which block/droplet a physical strand corresponds
  to) is tracked via a strand_id in the object's manifest, standing in
  for a real embedded barcode/index region within the oligo.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict, dataclass

from .addressing import PrimerRegistry
from .checksum import crc16_bytes
from .codec import get_codec, Strand, available_codecs
from .codec.naive import NaiveCodec
from .ecc.reed_solomon import ReedSolomonErasureCoder, ErasureBlock
from .formats import FastaRecord, FastqRecord, write_fasta, write_fastq, uniform_quality_string
from .metadata import MetadataStore
from .store import PhysicalPool
from .synthesis import SynthesisSimulator, SynthesisOrder, TwistSynthesisBackend
from .sequencing import SequencingSimulator, NanoporeSequencingBackend

_naive = NaiveCodec()
_CRC_BASE_LENGTH = len(_naive.encode(b"\x00\x00")[0])  # 8 bases for a 2-byte CRC


class DNAStorageError(Exception):
    pass


@dataclass
class VersionInfo:
    name: str
    version_number: int
    size_bytes: int
    codec: str
    ecc: str | None
    strand_count: int
    checksum: str


class DNAStorage:
    def __init__(
        self,
        backend: str = "simulator",
        codec: str = "fountain",
        parity_fraction: float = 0.4,
        fountain_redundancy: float = 2.5,
        state_path: str | None = None,
        pool_path: str | None = None,
        synthesis_kwargs: dict | None = None,
        sequencing_kwargs: dict | None = None,
    ):
        if codec not in available_codecs():
            raise ValueError(f"unknown codec '{codec}'. Available: {available_codecs()}")

        self.default_codec_name = codec
        self.parity_fraction = parity_fraction
        self.fountain_redundancy = fountain_redundancy
        self.state_path = state_path

        self.pool = PhysicalPool(pool_path)

        if backend == "simulator":
            self.synthesis = SynthesisSimulator(self.pool, **(synthesis_kwargs or {}))
            self.sequencer = SequencingSimulator(self.pool, **(sequencing_kwargs or {}))
        elif backend == "twist+nanopore":
            self.synthesis = TwistSynthesisBackend(**(synthesis_kwargs or {}))
            self.sequencer = NanoporeSequencingBackend(**(sequencing_kwargs or {}))
        else:
            raise ValueError(f"unknown backend '{backend}'")

        self.metadata = MetadataStore()
        self.primers = PrimerRegistry()
        self.manifests: dict[str, dict] = {}

        if state_path and os.path.exists(state_path):
            self._load_state()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def store(self, source: str | bytes, name: str | None = None, codec: str | None = None) -> VersionInfo:
        data, resolved_name = self._resolve_source(source, name)
        return self._write_version(resolved_name, data, codec)

    def update(self, source: str | bytes, name: str | None = None, codec: str | None = None) -> VersionInfo:
        """Semantically identical to store() -- included because DNA's
        append-only nature makes 'update' and 'create new version' the
        same physical operation. Both mark prior versions obsolete."""
        return self.store(source, name=name, codec=codec)

    def retrieve(self, name: str) -> bytes:
        record = self.metadata.get_by_name(name)
        if record is None:
            raise FileNotFoundError(f"no such object: '{name}'")
        version = record.latest_version()
        if version is None:
            raise FileNotFoundError(f"'{name}' has no active version")

        manifest = self.manifests[version.version_id]
        primer = self.primers.get(record.object_id)
        if primer is None:
            raise DNAStorageError(f"no primer pair registered for '{name}' -- archive state is inconsistent")

        read_results = self.sequencer.sequence(*primer)
        read_by_id = {r.strand_id: r.sequence for r in read_results}

        codec_obj = get_codec(version.codec)
        # (stripe_index, block_index) -> (payload_bases, is_parity)
        surviving: dict[tuple[int, int], tuple[str, bool]] = {}
        for rec in manifest["strand_records"]:
            seq = read_by_id.get(rec["strand_id"])
            if seq is None or len(seq) <= _CRC_BASE_LENGTH:
                continue  # dropped out, or too damaged by indels to even hold a CRC
            payload_bases, crc_bases = seq[:-_CRC_BASE_LENGTH], seq[-_CRC_BASE_LENGTH:]
            try:
                crc_recovered = _naive.decode([Strand(crc_bases)], 2)
                expected = crc16_bytes(payload_bases.encode("ascii"))
            except (KeyError, ValueError):
                continue  # unparseable -- treat as erased
            if crc_recovered != expected:
                continue  # corrupted -- treat as erased, let ECC/fountain recover it
            surviving[(rec["stripe_index"], rec["block_index"])] = (payload_bases, rec["is_parity"])

        if version.ecc == "reed_solomon":
            data = self._reconstruct_with_rs(name, manifest, codec_obj, surviving)
        else:
            strands = [Strand(payload) for (payload, _is_parity) in surviving.values()]
            data = codec_obj.decode(strands, manifest["original_length"])

        if hashlib.sha256(data).hexdigest() != version.checksum:
            raise DNAStorageError(
                f"checksum mismatch reconstructing '{name}': data recovered but does not "
                f"match the original checksum -- this should not happen if ECC/fountain "
                f"decode succeeded and indicates a bug rather than an uncorrected error"
            )
        return data

    def delete(self, name: str) -> None:
        """Tombstone the object. DNA is not physically erased -- use
        purge() for that."""
        self.metadata.tombstone(name)
        self._persist_state()

    def purge(self, name: str) -> int:
        """Physically discard all strand records for every version of
        `name`, and remove its metadata entirely. Returns the number of
        strand records removed. This models a real, deliberate 'destroy
        this sample' lab operation -- irreversible, unlike delete()."""
        record = self.metadata.purge(name)
        primer = self.primers.get(record.object_id)
        removed = 0
        if primer:
            removed = self.pool.purge_by_primer(*primer)
            self.primers.release(record.object_id)
        for version in record.versions:
            self.manifests.pop(version.version_id, None)
        self._persist_state()
        return removed

    def list_versions(self, name: str) -> list[VersionInfo]:
        record = self.metadata.get_by_name(name, include_tombstoned=True)
        if record is None:
            raise FileNotFoundError(f"no such object: '{name}'")
        return [
            VersionInfo(
                name=name, version_number=v.version_number, size_bytes=v.size_bytes,
                codec=v.codec, ecc=v.ecc, strand_count=v.strand_count, checksum=v.checksum,
            )
            for v in record.versions
        ]

    def list_objects(self, include_tombstoned: bool = False) -> list[str]:
        return [r.name for r in self.metadata.list_objects(include_tombstoned=include_tombstoned)]

    def stats(self) -> dict:
        return self.pool.stats()

    # ------------------------------------------------------------------
    # FASTA/FASTQ export (bioinformatics interop)
    # ------------------------------------------------------------------

    def export_object_fasta(self, name: str, path: str) -> int:
        """Write every strand belonging to `name`'s latest version as a
        multi-FASTA file -- the format real synthesis vendors actually
        accept for order submission. Returns the number of records
        written."""
        record = self.metadata.get_by_name(name)
        if record is None:
            raise FileNotFoundError(f"no such object: '{name}'")
        version = record.latest_version()
        if version is None:
            raise FileNotFoundError(f"'{name}' has no active version")
        manifest = self.manifests[version.version_id]
        primer = self.primers.get(record.object_id)

        fasta_records = []
        for rec in manifest["strand_records"]:
            sequence = self.pool.get(rec["strand_id"])
            if sequence is None:
                continue  # dropped at synthesis -- nothing to export for this strand
            desc = (
                f"object={name} stripe={rec['stripe_index']} block={rec['block_index']} "
                f"parity={rec['is_parity']} primer_fwd={primer[0]} primer_rev={primer[1]}"
            )
            fasta_records.append(FastaRecord(id=rec["strand_id"], sequence=sequence, description=desc))

        write_fasta(fasta_records, path)
        return len(fasta_records)

    def export_pool_fasta(self, path: str) -> int:
        """Write the entire physical pool (every strand for every object,
        including obsolete/tombstoned versions still physically present)
        as a single multi-FASTA file. Returns the number of records
        written."""
        fasta_records = [
            FastaRecord(id=strand_id, sequence=sequence, description="")
            for strand_id, sequence in self.pool.all_records()
            if sequence is not None
        ]
        write_fasta(fasta_records, path)
        return len(fasta_records)

    def retrieve_as_fastq(self, name: str, path: str, phred_score: int = 37) -> int:
        """Trigger a real sequencing read (through the configured
        SequencingBackend, with whatever error injection that implies)
        and dump the raw reads as FASTQ -- i.e. what a downstream
        bioinformatics pipeline would actually receive off a sequencer,
        before any DNAStorage-side error correction is applied. Returns
        the number of reads written.

        Quality scores are a flat placeholder (see
        formats.uniform_quality_string) since the simulator doesn't yet
        track per-base error provenance -- see CONTRIBUTING.md.
        """
        record = self.metadata.get_by_name(name)
        if record is None:
            raise FileNotFoundError(f"no such object: '{name}'")
        primer = self.primers.get(record.object_id)
        if primer is None:
            raise DNAStorageError(f"no primer pair registered for '{name}'")

        reads = self.sequencer.sequence(*primer)
        fastq_records = [
            FastqRecord(
                id=r.strand_id, sequence=r.sequence,
                quality=uniform_quality_string(len(r.sequence), phred_score),
                description=f"object={name}",
            )
            for r in reads
        ]
        write_fastq(fastq_records, path)
        return len(fastq_records)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _resolve_source(self, source: str | bytes, name: str | None) -> tuple[bytes, str]:
        if isinstance(source, (bytes, bytearray)):
            if name is None:
                raise ValueError("name= is required when storing raw bytes")
            return bytes(source), name
        if isinstance(source, str) and os.path.isfile(source):
            with open(source, "rb") as f:
                return f.read(), (name or os.path.basename(source))
        raise FileNotFoundError(f"'{source}' is not a file on disk; pass bytes with name= for in-memory data")

    # GF(256) caps a single Reed-Solomon codeword at 255 total shards.
    # Files needing more data shards than that are split into independent
    # RS-coded "stripes" (same principle as RAID striping across multiple
    # parity groups) so there's no hard file-size ceiling.
    STRIPE_MAX_DATA_SHARDS = 150

    def _write_version(self, name: str, data: bytes, codec_name: str | None) -> VersionInfo:
        codec_name = codec_name or self.default_codec_name
        # fountain_redundancy is tuning specific to our built-in fountain
        # codec's constructor; third-party codecs configure themselves via
        # their own defaults, so this kwarg is only passed for that one name.
        codec_obj = get_codec(codec_name, redundancy=self.fountain_redundancy) if codec_name == "fountain" else get_codec(codec_name)
        checksum_hex = hashlib.sha256(data).hexdigest()

        existing = self.metadata.get_by_name(name, include_tombstoned=True)
        if existing is None:
            existing = self.metadata.create_object(name)
        primer = self.primers.allocate(existing.object_id)

        version_number = len(existing.versions) + 1
        strand_records = []
        orders = []
        stripes_info = []

        if codec_obj.self_redundant:
            strands = codec_obj.encode(data)
            for i, strand in enumerate(strands):
                self._queue_strand(strand, existing.object_id, version_number, 0, i, False, primer, orders, strand_records)
            ecc_name = None
        else:
            capacity = codec_obj.shard_capacity_bytes()
            bytes_per_stripe = capacity * self.STRIPE_MAX_DATA_SHARDS
            num_stripes = max(1, -(-len(data) // bytes_per_stripe)) if data else 1
            ecc_name = "reed_solomon"

            for stripe_idx in range(num_stripes):
                stripe_data = data[stripe_idx * bytes_per_stripe:(stripe_idx + 1) * bytes_per_stripe]
                if not stripe_data and num_stripes > 1:
                    continue
                k = max(1, -(-len(stripe_data) // capacity)) if stripe_data else 1
                m = max(3, min(round(k * self.parity_fraction), 254 - k))
                rs = ReedSolomonErasureCoder(data_shards=k, parity_shards=m)
                blocks = rs.encode(stripe_data)
                stripes_info.append({"k": k, "m": m, "shard_len": len(blocks[0].payload), "length": len(stripe_data)})
                for b in blocks:
                    strand = codec_obj.encode(b.payload)[0]
                    self._queue_strand(
                        strand, existing.object_id, version_number, stripe_idx, b.index,
                        b.is_parity, primer, orders, strand_records,
                    )

        self.synthesis.synthesize(orders)

        version = self.metadata.add_version(
            name=name, size_bytes=len(data), checksum=checksum_hex, codec=codec_name,
            ecc=ecc_name, strand_count=len(strand_records),
        )
        self.manifests[version.version_id] = {
            "strand_records": strand_records,
            "stripes": stripes_info,  # empty for fountain
            "original_length": len(data),
        }
        self._persist_state()

        return VersionInfo(
            name=name, version_number=version.version_number, size_bytes=len(data),
            codec=codec_name, ecc=ecc_name, strand_count=len(strand_records), checksum=checksum_hex,
        )

    @staticmethod
    def _queue_strand(strand, object_id, version_number, stripe_index, block_index, is_parity, primer, orders, strand_records):
        crc_bases = _naive.encode(crc16_bytes(str(strand).encode("ascii")))[0]
        full_sequence = str(strand) + str(crc_bases)
        strand_id = f"{object_id}:{version_number}:{stripe_index}:{block_index}:{'p' if is_parity else 'd'}"
        orders.append(SynthesisOrder(strand_id, full_sequence, primer[0], primer[1]))
        strand_records.append({
            "strand_id": strand_id,
            "stripe_index": stripe_index,
            "block_index": block_index,
            "is_parity": is_parity,
        })

    def _reconstruct_with_rs(self, name, manifest, codec_obj, surviving) -> bytes:
        out = bytearray()
        for stripe_idx, ecc_info in enumerate(manifest["stripes"]):
            rs = ReedSolomonErasureCoder(data_shards=ecc_info["k"], parity_shards=ecc_info["m"])
            blocks = []
            for (s_idx, block_index), (payload_bases, is_parity) in surviving.items():
                if s_idx != stripe_idx:
                    continue
                payload_bytes = codec_obj.decode([Strand(payload_bases)], ecc_info["shard_len"])
                blocks.append(ErasureBlock(index=block_index, is_parity=is_parity, payload=payload_bytes))
            if len(blocks) < ecc_info["k"]:
                raise DNAStorageError(
                    f"insufficient surviving strands to reconstruct '{name}' (stripe {stripe_idx}): "
                    f"{len(blocks)}/{ecc_info['k']} needed -- too many strands were lost or corrupted"
                )
            out.extend(rs.decode(blocks, ecc_info["length"]))
        return bytes(out[:manifest["original_length"]])

    def _persist_state(self) -> None:
        if not self.state_path:
            return
        state = {
            "metadata": self.metadata.to_dict(),
            "primers": self.primers.to_dict(),
            "manifests": self.manifests,
        }
        with open(self.state_path, "w") as f:
            json.dump(state, f)

    def _load_state(self) -> None:
        with open(self.state_path, "r") as f:
            state = json.load(f)
        self.metadata = MetadataStore.from_dict(state.get("metadata", {}))
        self.primers = PrimerRegistry.from_dict(state.get("primers", {}))
        self.manifests = state.get("manifests", {})
