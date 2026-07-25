"""
Stub backend for a real sequencer (Oxford Nanopore's MinKNOW/API is
used as the example shape). Like the Twist synthesis stub, this is a
deliberately unimplemented integration seam: real sequencing runs
produce raw signal data that needs basecalling, demultiplexing by
primer/barcode, and quality filtering before you get back strand-like
reads -- a real implementation would sit on top of a basecalling
pipeline (e.g. Guppy/Dorado output), not just an HTTP client.
"""
from __future__ import annotations

import os

from .base import SequencingBackend, SequencedStrand


class NanoporeSequencingBackend(SequencingBackend):
    def __init__(self, device_address: str | None = None, api_token: str | None = None):
        self.device_address = device_address or os.environ.get("MINKNOW_DEVICE_ADDRESS")
        self.api_token = api_token or os.environ.get("MINKNOW_API_TOKEN")

    def is_available(self) -> bool:
        return bool(self.device_address)

    def sequence(self, primer_forward: str, primer_reverse: str) -> list[SequencedStrand]:
        if not self.is_available():
            raise RuntimeError(
                "NanoporeSequencingBackend requires a device address (set "
                "MINKNOW_DEVICE_ADDRESS or pass device_address=...). This backend "
                "is a stub: implement the MinKNOW API connection, run demultiplexing "
                "by primer_forward/primer_reverse against basecalled reads, and "
                "return SequencedStrand objects for matches before using it for "
                "real sequencing runs."
            )
        raise NotImplementedError(
            "Real sequencing readback is not implemented. This class exists as "
            "the integration point for a basecalling + demultiplexing pipeline."
        )
