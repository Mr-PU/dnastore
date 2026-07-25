"""
Synthesis simulator.

Models the error modes real DNA synthesis actually exhibits, so the
rest of the stack (codecs, ECC, addressing) can be developed and
benchmarked without lab access:

- substitution_rate: per-base chance a wrong nucleotide gets
  incorporated during synthesis.
- dropout_rate: per-strand chance synthesis fails outright for that
  oligo (this is the dominant real-world failure mode for DNA
  storage, dwarfing per-base substitution errors).

Indels during synthesis itself are rarer than during sequencing/PCR in
practice, so they're modeled in the sequencing simulator instead; see
`dnastore.sequencing.simulator`.
"""
from __future__ import annotations

import random

from .base import SynthesisBackend, SynthesisOrder
from ..codec.base import BASES


class SynthesisSimulator(SynthesisBackend):
    def __init__(self, pool, substitution_rate: float = 0.0002, dropout_rate: float = 0.01, seed: int | None = None):
        """`pool` is a dnastore.store.PhysicalPool instance shared with the
        matching SequencingSimulator."""
        self.pool = pool
        self.substitution_rate = substitution_rate
        self.dropout_rate = dropout_rate
        self._rng = random.Random(seed)

    def is_available(self) -> bool:
        return True

    def synthesize(self, orders: list[SynthesisOrder]) -> None:
        for order in orders:
            if self._rng.random() < self.dropout_rate:
                self.pool.write(order.strand_id, None, order.primer_forward, order.primer_reverse)
                continue

            sequence = self._apply_substitutions(order.sequence)
            self.pool.write(order.strand_id, sequence, order.primer_forward, order.primer_reverse)

    def _apply_substitutions(self, sequence: str) -> str:
        if self.substitution_rate <= 0:
            return sequence
        out = list(sequence)
        for i in range(len(out)):
            if self._rng.random() < self.substitution_rate:
                choices = [b for b in BASES if b != out[i]]
                out[i] = self._rng.choice(choices)
        return "".join(out)
