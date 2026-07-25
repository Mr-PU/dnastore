"""Small CRC16 (CCITT) implementation -- no external dependency needed.

Every DNA strand gets a checksum appended before synthesis. On
readback, a strand whose checksum doesn't match is treated as
*erased* rather than corrected in place -- correction is delegated to
the cross-strand Reed-Solomon erasure layer, which is a much better
fit for DNA's dominant failure modes (dropout, indels) than trying to
do single-strand error correction.
"""
from __future__ import annotations

_POLY = 0x1021
_INIT = 0xFFFF


def crc16(data: bytes) -> int:
    crc = _INIT
    for byte in data:
        crc ^= byte << 8
        for _ in range(8):
            if crc & 0x8000:
                crc = ((crc << 1) ^ _POLY) & 0xFFFF
            else:
                crc = (crc << 1) & 0xFFFF
    return crc


def crc16_bytes(data: bytes) -> bytes:
    return crc16(data).to_bytes(2, "big")


def verify(data: bytes, checksum: bytes) -> bool:
    return crc16_bytes(data) == checksum
