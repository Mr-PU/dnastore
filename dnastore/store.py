"""
PhysicalPool: the "physical DNA storage" layer.

In a real system this would be an actual tube/plate of synthesized DNA
molecules. Here it's the shared state that a SynthesisBackend writes
into and a SequencingBackend reads back out of, indexed by primer pair
so reads can be selective (random access) instead of "read
everything". Persisted to a JSON file so an archive survives across
process restarts, same as a real storage engine would persist its
on-disk representation.
"""
from __future__ import annotations

import json
import os
import threading


class PhysicalPool:
    def __init__(self, path: str | None = None):
        self.path = path
        self._lock = threading.Lock()
        # strand_id -> {sequence: str|None, primer_forward: str, primer_reverse: str}
        self._records: dict[str, dict] = {}
        if path and os.path.exists(path):
            self._load()

    def write(self, strand_id: str, sequence: str | None, primer_forward: str, primer_reverse: str) -> None:
        with self._lock:
            self._records[strand_id] = {
                "sequence": sequence,
                "primer_forward": primer_forward,
                "primer_reverse": primer_reverse,
            }
            self._maybe_persist()

    def read_by_primer(self, primer_forward: str, primer_reverse: str) -> list[tuple[str, str | None]]:
        with self._lock:
            return [
                (strand_id, rec["sequence"])
                for strand_id, rec in self._records.items()
                if rec["primer_forward"] == primer_forward and rec["primer_reverse"] == primer_reverse
            ]

    def get(self, strand_id: str) -> str | None:
        """Return the (possibly error-injected) sequence for a strand_id,
        or None if it was never written or dropped out at synthesis."""
        with self._lock:
            rec = self._records.get(strand_id)
            return rec["sequence"] if rec else None

    def all_records(self) -> list[tuple[str, str | None]]:
        """Return (strand_id, sequence) for every strand in the pool,
        across all objects/versions -- sequence is None for strands that
        dropped out at synthesis."""
        with self._lock:
            return [(strand_id, rec["sequence"]) for strand_id, rec in self._records.items()]

    def purge_by_primer(self, primer_forward: str, primer_reverse: str) -> int:
        """Physically discard all strands for a primer pair. Returns the
        number of strand records removed."""
        with self._lock:
            to_remove = [
                strand_id for strand_id, rec in self._records.items()
                if rec["primer_forward"] == primer_forward and rec["primer_reverse"] == primer_reverse
            ]
            for strand_id in to_remove:
                del self._records[strand_id]
            self._maybe_persist()
            return len(to_remove)

    def stats(self) -> dict:
        with self._lock:
            total = len(self._records)
            dropped = sum(1 for r in self._records.values() if r["sequence"] is None)
            return {"total_strands": total, "dropped_at_synthesis": dropped}

    def _maybe_persist(self) -> None:
        if self.path:
            self._save()

    def _save(self) -> None:
        with open(self.path, "w") as f:
            json.dump(self._records, f)

    def _load(self) -> None:
        with open(self.path, "r") as f:
            self._records = json.load(f)
