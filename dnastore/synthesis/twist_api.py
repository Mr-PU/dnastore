"""
Stub backend for a real DNA synthesis vendor (Twist Bioscience is used
as the example name/API shape; swap in whichever vendor you have lab
access through -- the interface is what matters).

This is intentionally unimplemented beyond request plumbing: real
synthesis is an async, multi-day turnaround process gated by account
approval and lab capacity, so a real integration needs a job-polling
or webhook flow this stub doesn't attempt to guess at. Treat this as
the seam where that integration would go, not a working client.
"""
from __future__ import annotations

import os

from .base import SynthesisBackend, SynthesisOrder


class TwistSynthesisBackend(SynthesisBackend):
    API_BASE_URL = "https://api.twistbioscience-vendor.example/v1"  # placeholder -- confirm with vendor docs

    def __init__(self, api_key: str | None = None):
        self.api_key = api_key or os.environ.get("TWIST_API_KEY")

    def is_available(self) -> bool:
        return bool(self.api_key)

    def synthesize(self, orders: list[SynthesisOrder]) -> None:
        if not self.is_available():
            raise RuntimeError(
                "TwistSynthesisBackend requires an API key (set TWIST_API_KEY or pass "
                "api_key=...). This backend is a stub: implement submit_order()/"
                "poll_job() against your vendor's actual API before using it for "
                "real synthesis orders."
            )
        raise NotImplementedError(
            "Real synthesis submission is not implemented. This class exists as "
            "the integration point: implement HTTP calls to your vendor's order "
            "and job-status endpoints here, store the returned job ID per strand, "
            "and poll until synthesis completes before considering strands written."
        )
