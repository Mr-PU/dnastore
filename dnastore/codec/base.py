"""Codec interface: binary bytes <-> DNA base strings.

Every codec turns a bytes payload into one or more DNA "strands"
(each a string over {A, C, G, T} bounded by MAX_STRAND_LENGTH, which
mirrors realistic synthesis constraints) and back. Codecs are
pluggable so naive / rotating / fountain encodings can be benchmarked
against the same error injector.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

BASES = "ACGT"
BASE_TO_BITS = {"A": 0b00, "C": 0b01, "G": 0b10, "T": 0b11}
BITS_TO_BASE = {v: k for k, v in BASE_TO_BITS.items()}

# Real synthesis vendors generally top out well under 300 nt per oligo.
MAX_STRAND_LENGTH = 200


class Strand(str):
    """A DNA strand is just a string over {A,C,G,T}; this subclass exists
    purely for type clarity in signatures."""


class Codec(ABC):
    name: str = "base"

    # Whether this codec supplies its own cross-strand redundancy (e.g. a
    # rateless fountain code) and should therefore skip DNAStorage's
    # external Reed-Solomon erasure layer. False (the default) means
    # "wrap this codec's output in RS erasure coding" -- override to True
    # for self-redundant codecs.
    self_redundant: bool = False

    # Optional explicit override for shard_capacity_bytes(). Concrete
    # codecs can set this class attribute for speed; if left None, the
    # default shard_capacity_bytes() implementation auto-detects it by
    # probing encode(), which works for any well-behaved codec (including
    # third-party plugins) without requiring them to know this exists.
    _shard_capacity_bytes: int | None = None

    @abstractmethod
    def encode(self, data: bytes) -> list[Strand]:
        """Encode a bytes payload into one or more DNA strands."""

    @abstractmethod
    def decode(self, strands: list[Strand], original_length: int) -> bytes:
        """Decode DNA strands back into the original bytes payload.

        `strands` may be a subset/superset depending on codec (fountain
        codecs tolerate missing/extra droplets; naive/rotating codecs
        expect all strands present, in order, with any redundancy
        handled by the ECC layer instead).
        """

    def shard_capacity_bytes(self) -> int:
        """The largest number of input bytes this codec is guaranteed to
        pack into exactly one DNA strand. DNAStorage uses this to size
        Reed-Solomon shards for non-self-redundant codecs, so that each
        RS shard maps 1:1 onto one physical strand.

        Auto-detected by probing encode() if a codec doesn't override
        `_shard_capacity_bytes` -- third-party plugin codecs don't need
        to know this mechanism exists for DNAStorage to work correctly,
        though declaring `_shard_capacity_bytes` explicitly is faster
        (skips the probe) and recommended once you know the number.
        """
        if self._shard_capacity_bytes is not None:
            return self._shard_capacity_bytes

        n = 1
        last_good = 1
        while n < 100_000 and len(self.encode(bytes(n))) <= 1:
            last_good = n
            n *= 2
        lo, hi = last_good, n
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if len(self.encode(bytes(mid))) <= 1:
                lo = mid
            else:
                hi = mid - 1
        return max(1, lo)

    @staticmethod
    def gc_content(strand: str) -> float:
        if not strand:
            return 0.0
        gc = sum(1 for b in strand if b in "GC")
        return gc / len(strand)

    @staticmethod
    def max_homopolymer_run(strand: str) -> int:
        if not strand:
            return 0
        best = cur = 1
        for i in range(1, len(strand)):
            if strand[i] == strand[i - 1]:
                cur += 1
                best = max(best, cur)
            else:
                cur = 1
        return best
