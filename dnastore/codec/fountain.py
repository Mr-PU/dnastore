"""
DNA Fountain codec, following the rateless-code approach popularised
for DNA storage by Erlich & Zielinski (2017, "DNA Fountain enables a
robust and efficient storage architecture").

Unlike naive/rotating (which rely on the external Reed-Solomon
erasure layer for redundancy), a fountain code is rateless: the
encoder can produce an unlimited stream of "droplets", each formed by
XOR-ing a random subset of the K source blocks. The subset for a
droplet is *not stored explicitly* -- it's regenerated from a small
integer seed using a seeded PRNG, so the only overhead per droplet is
that seed. The decoder regenerates the same neighbor sets from the
same seeds and reconstructs the K source blocks via belief-propagation
("peeling"): resolve any droplet that currently has exactly one
unresolved neighbor, subtract it out, and repeat.

This makes fountain codes naturally tolerant of strand dropout: as
long as *some* sufficiently diverse (K * ~1.05-1.5, depending on
degree distribution quality) droplets survive -- regardless of which
ones -- decoding succeeds. No separate erasure-coding layer needed.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass

from .base import Codec, Strand
from .naive import NaiveCodec

_naive = NaiveCodec()

SEED_BYTES = 4
DEFAULT_REDUNDANCY = 1.6  # produce ~60% more droplets than source blocks


class FountainDecodeError(Exception):
    pass


@dataclass
class _Droplet:
    seed: int
    payload: bytes


_soliton_cache: dict[tuple, list[float]] = {}


def _robust_soliton_cdf(k: int, c: float = 0.2, delta: float = 0.05) -> list[float]:
    """Cumulative distribution (length k, index d-1 -> P(degree <= d)) for
    Luby's robust soliton distribution. This adds a small spike around
    degree k/S on top of the ideal soliton distribution rho(), which is
    what actually guarantees (with high probability) that enough
    degree-1 droplets exist throughout peeling -- the ideal soliton
    alone fails far too often in practice, which is why real fountain-
    code implementations (including DNA Fountain) use the robust
    version instead.
    """
    key = (k, c, delta)
    if key in _soliton_cache:
        return _soliton_cache[key]

    rho = [0.0] * (k + 1)  # 1-indexed
    rho[1] = 1.0 / k
    for d in range(2, k + 1):
        rho[d] = 1.0 / (d * (d - 1))

    S = max(1.0, c * math.sqrt(k) * math.log(max(k, 2) / delta))
    spike = max(1, min(k - 1, round(k / S)))

    tau = [0.0] * (k + 1)
    for d in range(1, spike):
        tau[d] = S / (k * d)
    tau[spike] = S * math.log(S / delta) / k

    combined = [rho[d] + tau[d] for d in range(k + 1)]
    total = sum(combined[1:])
    mu = [0.0] + [v / total for v in combined[1:]]

    cdf = []
    running = 0.0
    for d in range(1, k + 1):
        running += mu[d]
        cdf.append(running)
    cdf[-1] = 1.0  # guard against float drift
    _soliton_cache[key] = cdf
    return cdf


def _sample_degree(cdf: list[float], rng: random.Random) -> int:
    u = rng.random()
    lo, hi = 0, len(cdf) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if cdf[mid] < u:
            lo = mid + 1
        else:
            hi = mid
    return lo + 1


class FountainCodec(Codec):
    name = "fountain"
    self_redundant = True

    def __init__(self, segment_bytes: int = 32, redundancy: float = DEFAULT_REDUNDANCY):
        self.segment_bytes = segment_bytes
        self.redundancy = redundancy

    def encode(self, data: bytes) -> list[Strand]:
        k = max(1, -(-len(data) // self.segment_bytes))
        padded = data + b"\x00" * (k * self.segment_bytes - len(data))
        blocks = [padded[i * self.segment_bytes:(i + 1) * self.segment_bytes] for i in range(k)]

        cdf = _robust_soliton_cdf(k)

        def make_droplet() -> tuple[int, list[int], bytearray]:
            seed = random.randrange(2**32)
            rng = random.Random(seed)
            degree = _sample_degree(cdf, rng)
            neighbors = rng.sample(range(k), degree)
            payload = bytearray(self.segment_bytes)
            for idx in neighbors:
                block = blocks[idx]
                for i in range(self.segment_bytes):
                    payload[i] ^= block[i]
            return seed, neighbors, payload

        droplets: list[tuple[int, list[int], bytearray]] = []
        batch_size = max(1, int(k * self.redundancy) - k) or max(1, k // 4)
        target = max(k, int(k * self.redundancy))

        # Generate an initial batch, then keep growing until a self-check
        # peeling decode (no simulated transmission errors -- this only
        # verifies the *code* has enough degree diversity, not the
        # physical channel) confirms every block is recoverable. This
        # trades a little extra compute at encode time for a hard
        # guarantee instead of hoping a fixed multiplier was enough.
        max_droplets = k * 6 + 50
        while True:
            while len(droplets) < target:
                droplets.append(make_droplet())
            if self._peelable(k, [(s, n, p) for s, n, p in droplets]):
                break
            if len(droplets) >= max_droplets:
                raise FountainDecodeError(
                    f"could not construct a recoverable fountain code for k={k} "
                    f"blocks even after {len(droplets)} droplets; try a larger "
                    f"segment_bytes or a different redundancy"
                )
            target = len(droplets) + batch_size

        strands: list[Strand] = []
        for seed, _neighbors, payload in droplets:
            header_bases = _naive.encode(seed.to_bytes(SEED_BYTES, "big"))[0]
            payload_bases = _naive.encode(bytes(payload))[0]
            strands.append(Strand(str(header_bases) + str(payload_bases)))
        return strands

    @staticmethod
    def _peelable(k: int, droplets: list[tuple[int, list[int], bytearray]]) -> bool:
        """Cheap self-check: can these droplets' neighbor sets alone (index
        sets, ignoring actual payload correctness which is guaranteed by
        construction) be peeled down to all k blocks?"""
        pending = [set(neighbors) for _, neighbors, _ in droplets]
        resolved: set[int] = set()
        index_to_entries: dict[int, list[int]] = {i: [] for i in range(k)}
        for pos, neighbors in enumerate(pending):
            for idx in neighbors:
                index_to_entries[idx].append(pos)

        queue = [pos for pos, n in enumerate(pending) if len(n) == 1]
        while queue and len(resolved) < k:
            pos = queue.pop()
            neighbors = pending[pos]
            if len(neighbors) != 1:
                continue
            idx = next(iter(neighbors))
            if idx in resolved:
                continue
            resolved.add(idx)
            neighbors.clear()
            for other_pos in index_to_entries[idx]:
                other = pending[other_pos]
                if idx in other:
                    other.discard(idx)
                    if len(other) == 1:
                        queue.append(other_pos)
        return len(resolved) >= k

    def decode(self, strands: list[Strand], original_length: int) -> bytes:
        k = max(1, -(-original_length // self.segment_bytes))
        header_len_bases = len(_naive.encode(bytes(SEED_BYTES))[0])

        droplets: list[_Droplet] = []
        for strand in strands:
            header = strand[:header_len_bases]
            payload_bases = strand[header_len_bases:]
            seed_bytes = _naive.decode([Strand(header)], SEED_BYTES)
            seed = int.from_bytes(seed_bytes, "big")
            payload = _naive.decode([Strand(payload_bases)], self.segment_bytes)
            droplets.append(_Droplet(seed=seed, payload=payload))

        resolved: dict[int, bytearray] = {}

        # Each pending entry: [set of unresolved neighbor indices, current
        # XOR payload]. Invariant: payload always reflects the XOR of only
        # the blocks still in `neighbors` -- any resolved block is XORed
        # out the instant it's discovered, never left stale.
        cdf = _robust_soliton_cdf(k)
        pending: list[list] = []
        for d in droplets:
            rng = random.Random(d.seed)
            degree = _sample_degree(cdf, rng)
            neighbors = set(rng.sample(range(k), degree))
            pending.append([neighbors, bytearray(d.payload)])

        # index -> positions of pending entries that reference it, so
        # resolving a block only touches droplets that contain it.
        index_to_entries: dict[int, list[int]] = {i: [] for i in range(k)}
        for pos, (neighbors, _payload) in enumerate(pending):
            for idx in neighbors:
                index_to_entries[idx].append(pos)

        def _apply_resolution(idx: int, value: bytearray) -> list[int]:
            newly_ready = []
            for pos in index_to_entries[idx]:
                neighbors, payload = pending[pos]
                if idx not in neighbors:
                    continue
                for i in range(len(payload)):
                    payload[i] ^= value[i]
                neighbors.discard(idx)
                if len(neighbors) == 1:
                    newly_ready.append(pos)
            return newly_ready

        queue = [pos for pos, (neighbors, _p) in enumerate(pending) if len(neighbors) == 1]
        while queue and len(resolved) < k:
            pos = queue.pop()
            neighbors, payload = pending[pos]
            if len(neighbors) != 1:
                continue
            idx = next(iter(neighbors))
            if idx in resolved:
                continue
            resolved[idx] = bytearray(payload)
            neighbors.clear()
            queue.extend(_apply_resolution(idx, resolved[idx]))

        if len(resolved) < k:
            raise FountainDecodeError(
                f"decode failed: only resolved {len(resolved)}/{k} source blocks "
                f"from {len(strands)} droplets -- too many strands were lost/corrupted"
            )

        out = bytearray()
        for i in range(k):
            out.extend(resolved[i])
        return bytes(out[:original_length])
