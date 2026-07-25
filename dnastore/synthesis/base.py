"""Synthesis backend interface: turns DNA base-strings into physically
(or virtually) stored strand records."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class SynthesisOrder:
    strand_id: str
    sequence: str
    primer_forward: str
    primer_reverse: str


class SynthesisBackend(ABC):
    @abstractmethod
    def synthesize(self, orders: list[SynthesisOrder]) -> None:
        """Submit strands for synthesis. For the simulator this is
        effectively instantaneous and stores strands in-memory/on-disk
        with injected errors; for a real vendor backend this would be an
        async job (synthesis takes days), so a production implementation
        of a real backend would need polling/webhook support layered on
        top of this call.
        """

    @abstractmethod
    def is_available(self) -> bool:
        """Whether this backend is actually usable right now (e.g. has
        valid API credentials configured)."""
