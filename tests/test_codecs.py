import os
import pytest

from dnastore.codec import NaiveCodec, RotatingCodec, FountainCodec, get_codec


CODECS = [NaiveCodec(), RotatingCodec(), FountainCodec()]


@pytest.mark.parametrize("codec", CODECS, ids=lambda c: c.name)
def test_round_trip_random_bytes(codec):
    data = os.urandom(500)
    strands = codec.encode(data)
    decoded = codec.decode(strands, len(data))
    assert decoded == data


@pytest.mark.parametrize("codec", CODECS, ids=lambda c: c.name)
def test_round_trip_empty(codec):
    data = b""
    strands = codec.encode(data)
    decoded = codec.decode(strands, len(data))
    assert decoded == data


@pytest.mark.parametrize("codec", CODECS, ids=lambda c: c.name)
def test_round_trip_small(codec):
    data = b"A"
    strands = codec.encode(data)
    decoded = codec.decode(strands, len(data))
    assert decoded == data


def test_rotating_codec_has_no_homopolymers():
    codec = RotatingCodec()
    data = bytes([0] * 1000)  # worst case for naive encoding
    strands = codec.encode(data)
    for strand in strands:
        assert codec.max_homopolymer_run(strand) == 1


def test_naive_codec_can_produce_homopolymers():
    codec = NaiveCodec()
    data = bytes([0] * 100)  # all-zero bytes -> "AAAA..." under naive mapping
    strands = codec.encode(data)
    assert any(codec.max_homopolymer_run(s) > 1 for s in strands)


def test_fountain_tolerates_dropout():
    import random
    codec = FountainCodec(redundancy=2.5)
    data = os.urandom(2000)
    strands = codec.encode(data)
    surviving = [s for s in strands if random.random() > 0.1]
    decoded = codec.decode(surviving, len(data))
    assert decoded == data


def test_get_codec_registry():
    for name in ("naive", "rotating", "fountain"):
        codec = get_codec(name)
        assert codec.name == name


def test_get_codec_unknown_raises():
    with pytest.raises(ValueError):
        get_codec("not-a-real-codec")
