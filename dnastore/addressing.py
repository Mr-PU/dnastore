"""
Primer-based addressing.

Real DNA storage systems don't sequence the entire pool to read one
file back -- they flank each object's strands with a unique PCR
primer pair, then selectively amplify just that primer's strands
before sequencing. This module simulates that addressing scheme: it
allocates a unique primer pair per stored object and is what lets
`retrieve()` be a targeted read instead of "decode everything".

Primers are modeled as short DNA sequences here; a real implementation
would need to check new primers against melting temperature and
cross-hybridization constraints versus every existing primer in the
pool, which we approximate by simply guaranteeing uniqueness and a
minimum length.
"""
from __future__ import annotations

import hashlib

PRIMER_LENGTH = 20
PRIMER_BASES = "ACGT"


def _hash_to_primer(seed: str, length: int = PRIMER_LENGTH) -> str:
    digest = hashlib.sha256(seed.encode("utf-8")).digest()
    bits = "".join(f"{byte:08b}" for byte in digest)
    bases = []
    for i in range(0, length * 2, 2):
        pair = bits[i:i + 2]
        bases.append(PRIMER_BASES[int(pair, 2)])
    return "".join(bases)


class PrimerRegistry:
    """Allocates unique forward/reverse primer pairs per object."""

    def __init__(self):
        self._used_primers: set[str] = set()
        self._object_primers: dict[str, tuple[str, str]] = {}

    def allocate(self, object_id: str) -> tuple[str, str]:
        if object_id in self._object_primers:
            return self._object_primers[object_id]

        attempt = 0
        while True:
            forward = _hash_to_primer(f"{object_id}:fwd:{attempt}")
            reverse = _hash_to_primer(f"{object_id}:rev:{attempt}")
            if forward not in self._used_primers and reverse not in self._used_primers and forward != reverse:
                break
            attempt += 1

        self._used_primers.add(forward)
        self._used_primers.add(reverse)
        self._object_primers[object_id] = (forward, reverse)
        return forward, reverse

    def get(self, object_id: str) -> tuple[str, str] | None:
        return self._object_primers.get(object_id)

    def release(self, object_id: str) -> None:
        """Called on purge() -- frees the primer pair for potential reuse.
        Not called on delete()/tombstone, since a tombstoned object may
        still need its primers for audit/undelete purposes."""
        pair = self._object_primers.pop(object_id, None)
        if pair:
            self._used_primers.discard(pair[0])
            self._used_primers.discard(pair[1])

    def to_dict(self) -> dict:
        return {oid: list(pair) for oid, pair in self._object_primers.items()}

    @classmethod
    def from_dict(cls, data: dict) -> "PrimerRegistry":
        reg = cls()
        for oid, pair in data.items():
            reg._object_primers[oid] = tuple(pair)
            reg._used_primers.update(pair)
        return reg
