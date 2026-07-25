"""
Rotating ternary codec, in the spirit of Goldman et al. (2013).

Naive 2-bit encoding can produce long homopolymer runs (e.g. a run of
zero bytes -> "AAAA..."), which both synthesis and sequencing handle
badly. The fix used by early DNA-storage papers: encode data in base
3 (trits) instead of base 2, and map each trit to a DNA base chosen
from the *3 letters other than the previous base*. Structurally, this
makes a homopolymer run impossible -- every base differs from its
predecessor by construction, at the cost of some information density
(a trit carries log2(3) ~= 1.585 bits instead of 2).

Data is processed in fixed-size byte chunks, each converted to a
fixed number of trits (sized so 3^n_trits comfortably covers
256^chunk_bytes), then mapped through the rotation table.
"""
from __future__ import annotations

import math

from .base import BASES, Codec, Strand, MAX_STRAND_LENGTH

CHUNK_BYTES = 32
TRITS_PER_CHUNK = math.ceil(CHUNK_BYTES * 8 * math.log(2) / math.log(3))
assert TRITS_PER_CHUNK <= MAX_STRAND_LENGTH, "chunk too large for MAX_STRAND_LENGTH"

# Fixed pseudo-previous-base used at the start of every strand, so
# encode/decode agree without needing to share extra state.
_START_CONTEXT = "T"


def _rotation_table() -> dict[str, str]:
    """For each possible previous base, the 3 remaining bases in a fixed
    cyclic order (trit value 0/1/2 -> letter)."""
    table = {}
    for prev in BASES:
        table[prev] = [b for b in BASES if b != prev]
    return table


_ROTATE = _rotation_table()
_REVERSE_ROTATE = {prev: {letter: i for i, letter in enumerate(letters)} for prev, letters in _ROTATE.items()}


class RotatingCodec(Codec):
    name = "rotating"
    _shard_capacity_bytes = CHUNK_BYTES

    def encode(self, data: bytes) -> list[Strand]:
        strands: list[Strand] = []
        for start in range(0, max(len(data), 1), CHUNK_BYTES):
            chunk = data[start:start + CHUNK_BYTES]
            if not chunk and strands:
                break
            padded = chunk + b"\x00" * (CHUNK_BYTES - len(chunk))
            trits = self._chunk_to_trits(padded)
            strands.append(Strand(self._trits_to_strand(trits)))
            if not data:
                break
        return strands or [Strand(self._trits_to_strand([0] * TRITS_PER_CHUNK))]

    def decode(self, strands: list[Strand], original_length: int) -> bytes:
        out = bytearray()
        for strand in strands:
            trits = self._strand_to_trits(strand)
            out.extend(self._trits_to_chunk(trits))
        return bytes(out[:original_length])

    @staticmethod
    def _chunk_to_trits(chunk: bytes) -> list[int]:
        n = int.from_bytes(chunk, "big")
        trits = []
        for _ in range(TRITS_PER_CHUNK):
            trits.append(n % 3)
            n //= 3
        trits.reverse()
        return trits

    @staticmethod
    def _trits_to_chunk(trits: list[int]) -> bytes:
        n = 0
        for t in trits:
            n = n * 3 + t
        return n.to_bytes(CHUNK_BYTES, "big")

    @staticmethod
    def _trits_to_strand(trits: list[int]) -> str:
        out = []
        prev = _START_CONTEXT
        for t in trits:
            letter = _ROTATE[prev][t]
            out.append(letter)
            prev = letter
        return "".join(out)

    @staticmethod
    def _strand_to_trits(strand: str) -> list[int]:
        trits = []
        prev = _START_CONTEXT
        for letter in strand:
            trits.append(_REVERSE_ROTATE[prev][letter])
            prev = letter
        return trits
