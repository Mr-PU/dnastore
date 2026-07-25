"""Naive 2-bit-per-base codec: A=00, C=01, G=10, T=11.

This is the baseline every DNA-storage paper compares against. It's
maximally information-dense (2 bits/base) but produces long
homopolymer runs (e.g. all-zero bytes -> "AAAA...") and arbitrary
GC-content, both of which synthesis and sequencing handle poorly in
practice. Included so benchmarks have a naive control to beat.
"""
from __future__ import annotations

from .base import BASES, BITS_TO_BASE, BASE_TO_BITS, Codec, Strand, MAX_STRAND_LENGTH

BITS_PER_STRAND = MAX_STRAND_LENGTH * 2
BYTES_PER_STRAND = BITS_PER_STRAND // 8


class NaiveCodec(Codec):
    name = "naive"
    _shard_capacity_bytes = BYTES_PER_STRAND

    def encode(self, data: bytes) -> list[Strand]:
        strands: list[Strand] = []
        for chunk_start in range(0, len(data), BYTES_PER_STRAND):
            chunk = data[chunk_start:chunk_start + BYTES_PER_STRAND]
            strands.append(Strand(self._bytes_to_bases(chunk)))
        return strands or [Strand("")]

    def decode(self, strands: list[Strand], original_length: int) -> bytes:
        out = bytearray()
        for strand in strands:
            out.extend(self._bases_to_bytes(strand))
        return bytes(out[:original_length])

    @staticmethod
    def _bytes_to_bases(chunk: bytes) -> str:
        bases = []
        for byte in chunk:
            for shift in (6, 4, 2, 0):
                bases.append(BITS_TO_BASE[(byte >> shift) & 0b11])
        return "".join(bases)

    @staticmethod
    def _bases_to_bytes(strand: str) -> bytes:
        out = bytearray()
        bits = [BASE_TO_BITS[b] for b in strand]
        for i in range(0, len(bits) - len(bits) % 4, 4):
            byte = 0
            for j in range(4):
                byte = (byte << 2) | bits[i + j]
            out.append(byte)
        return bytes(out)
