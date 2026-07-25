"""
Reed-Solomon erasure coding across DNA strands.

DNA storage's dominant failure mode is not per-symbol substitution
inside a strand (the codec layer + per-strand checksum handles that by
treating a corrupted strand as "erased") -- it's whole-strand dropout:
some fraction of strands simply never come back from synthesis,
long-term storage, or PCR amplification.

This module implements systematic Reed-Solomon *erasure* coding over
GF(256): given k data strands, it produces m additional parity
strands. Any k of the resulting (k + m) strands are sufficient to
reconstruct the original data, regardless of *which* strands were
lost. This is the same principle used by real-world archival storage
(e.g. RAID 6, distributed object stores) and by DNA-storage papers
that use RS across oligos for exactly this reason.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..gf256 import Matrix, mul


@dataclass
class ErasureBlock:
    index: int
    is_parity: bool
    payload: bytes


class ReedSolomonErasureCoder:
    """
    Uses a full (not systematic) Vandermonde generator matrix: row i is
    [x_i^0, x_i^1, ..., x_i^(k-1)] for distinct nonzero x_i = i+1. This
    is deliberate, not a simplification -- it's what actually guarantees
    the MDS (maximum distance separable) property that *any* k of the n
    rows form an invertible k x k submatrix, because a square Vandermonde
    submatrix built from distinct evaluation points always has a nonzero
    determinant (the classic Vandermonde determinant formula). A
    systematic construction (identity rows for data + a separately-built
    Vandermonde tail for parity) does *not* automatically have this
    property -- mixing identity rows with unrelated Vandermonde rows can
    produce a singular submatrix for some erasure patterns, which is a
    real bug this implementation specifically avoids.
    """

    def __init__(self, data_shards: int, parity_shards: int):
        if data_shards <= 0 or parity_shards < 0:
            raise ValueError("data_shards must be > 0 and parity_shards >= 0")
        if data_shards + parity_shards > 255:
            raise ValueError("GF(256) supports at most 255 total shards")
        self.k = data_shards
        self.m = parity_shards
        self.n = data_shards + parity_shards
        self._generator = Matrix.vandermonde(self.n, self.k)

    def encode(self, data: bytes) -> list[ErasureBlock]:
        """Split `data` into k equal-length shards (zero-padded) and
        produce n = k + m encoded blocks via the Vandermonde generator.
        Blocks are not raw copies of the input (non-systematic), so
        decode() is required to recover the original data even when
        all k "data-labeled" blocks happen to survive."""
        shard_len = -(-len(data) // self.k) or 1  # ceil div, min 1
        padded = data + b"\x00" * (shard_len * self.k - len(data))
        shards = [padded[i * shard_len:(i + 1) * shard_len] for i in range(self.k)]

        blocks: list[ErasureBlock] = []
        for row in range(self.n):
            gen_row = self._generator.rows[row]
            encoded = bytearray(shard_len)
            for byte_pos in range(shard_len):
                acc = 0
                for j in range(self.k):
                    acc ^= mul(gen_row[j], shards[j][byte_pos])
                encoded[byte_pos] = acc
            blocks.append(ErasureBlock(index=row, is_parity=(row >= self.k), payload=bytes(encoded)))

        return blocks

    def decode(self, surviving: list[ErasureBlock], original_length: int) -> bytes:
        """Reconstruct the original bytes from any k surviving blocks."""
        if len(surviving) < self.k:
            raise ValueError(
                f"cannot reconstruct: need >= {self.k} surviving strands, got {len(surviving)}"
            )
        surviving = sorted(surviving, key=lambda b: b.index)[: self.k]
        shard_len = len(surviving[0].payload)

        used_rows = [b.index for b in surviving]
        sub_gen = self._generator.sub_matrix(used_rows)
        inv_gen = sub_gen.invert()

        received = Matrix([[b.payload[byte_pos] for byte_pos in range(shard_len)] for b in surviving])
        recovered = inv_gen.multiply(received)

        out = bytearray()
        for row in recovered.rows:
            out.extend(row)
        return bytes(out[:original_length])
