"""
A deliberately trivial codec, to demonstrate the plugin mechanism --
not a serious encoding scheme. It XORs every byte with a fixed key
before delegating to dnastore's built-in NaiveCodec for the actual
binary-to-DNA mapping.

This file, and the tiny pyproject.toml next to it, are the complete
template for a third-party dnastore codec plugin: implement
dnastore.codec.base.Codec, declare it under the "dnastore.codecs"
entry-point group, `pip install` your package, and it shows up in
dnastore.available_codecs() with zero changes to the dnastore repo.
"""
from __future__ import annotations

from dnastore.codec.base import Codec, Strand
from dnastore.codec.naive import NaiveCodec

_XOR_KEY = 0x5A
_naive = NaiveCodec()


class XorDemoCodec(Codec):
    name = "xor_demo"

    def encode(self, data: bytes) -> list[Strand]:
        transformed = bytes(b ^ _XOR_KEY for b in data)
        return _naive.encode(transformed)

    def decode(self, strands: list[Strand], original_length: int) -> bytes:
        transformed = _naive.decode(strands, original_length)
        return bytes(b ^ _XOR_KEY for b in transformed)
