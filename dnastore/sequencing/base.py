"""Sequencing backend interface: reads strands back out of physical
(or simulated) DNA storage, addressed by primer pair."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class SequencedStrand:
    strand_id: str
    sequence: str  # as read -- may contain errors relative to what was synthesized


class SequencingBackend(ABC):
    @abstractmethod
    def sequence(self, primer_forward: str, primer_reverse: str) -> list[SequencedStrand]:
        """Selectively amplify + read back strands matching a primer pair.
        Strands that dropped out during synthesis, or that fail to
        amplify/read here, are simply absent from the returned list --
        callers must be able to reconstruct from a subset (that's what
        the ECC/fountain layers are for)."""

    @abstractmethod
    def is_available(self) -> bool:
        ...
