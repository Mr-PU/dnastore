"""
Sequencing simulator.

Reads strands back from a shared PhysicalPool (populated by a
SynthesisSimulator), applying read-time error modes:

- substitution_rate / insertion_rate / deletion_rate: per-base error
  probabilities during sequencing/PCR. Indels matter more here than at
  synthesis time -- this is deliberately where the simulator models
  them, since real sequencing platforms (especially long-read/
  nanopore-style) are indel-prone in a way short-read synthesis
  generally isn't.
- read_dropout_rate: chance a strand that *was* successfully
  synthesized still fails to be read back (degradation, failed PCR
  amplification, sampling loss from a finite pool draw).
"""
from __future__ import annotations

import random

from .base import SequencingBackend, SequencedStrand
from ..codec.base import BASES


class SequencingSimulator(SequencingBackend):
    def __init__(
        self, pool, substitution_rate: float = 0.0002, insertion_rate: float = 0.00005,
        deletion_rate: float = 0.00005, read_dropout_rate: float = 0.01, seed: int | None = None,
    ):
        self.pool = pool
        self.substitution_rate = substitution_rate
        self.insertion_rate = insertion_rate
        self.deletion_rate = deletion_rate
        self.read_dropout_rate = read_dropout_rate
        self._rng = random.Random(seed)

    def is_available(self) -> bool:
        return True

    def sequence(self, primer_forward: str, primer_reverse: str) -> list[SequencedStrand]:
        records = self.pool.read_by_primer(primer_forward, primer_reverse)
        results: list[SequencedStrand] = []
        for strand_id, sequence in records:
            if sequence is None:  # dropped out at synthesis time
                continue
            if self._rng.random() < self.read_dropout_rate:
                continue
            read = self._apply_read_errors(sequence)
            results.append(SequencedStrand(strand_id=strand_id, sequence=read))
        return results

    def _apply_read_errors(self, sequence: str) -> str:
        out = []
        for base in sequence:
            r = self._rng.random()
            if r < self.deletion_rate:
                continue  # base dropped
            r2 = self._rng.random()
            if r2 < self.insertion_rate:
                out.append(self._rng.choice(BASES))  # extra inserted base before this one
            if self._rng.random() < self.substitution_rate:
                out.append(self._rng.choice([b for b in BASES if b != base]))
            else:
                out.append(base)
        return "".join(out)
