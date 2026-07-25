import os
import random
import pytest

from dnastore.ecc.reed_solomon import ReedSolomonErasureCoder


def test_encode_decode_no_loss():
    data = os.urandom(1000)
    rs = ReedSolomonErasureCoder(data_shards=10, parity_shards=4)
    blocks = rs.encode(data)
    assert len(blocks) == 14
    decoded = rs.decode(blocks, len(data))
    assert decoded == data


def test_decode_tolerates_max_erasures():
    data = os.urandom(1000)
    rs = ReedSolomonErasureCoder(data_shards=10, parity_shards=4)
    blocks = rs.encode(data)
    random.shuffle(blocks)
    surviving = blocks[4:]  # drop exactly `parity_shards` blocks
    decoded = rs.decode(surviving, len(data))
    assert decoded == data


def test_decode_fails_below_threshold():
    data = os.urandom(1000)
    rs = ReedSolomonErasureCoder(data_shards=10, parity_shards=4)
    blocks = rs.encode(data)
    random.shuffle(blocks)
    surviving = blocks[5:]  # drop one more than tolerable
    with pytest.raises(ValueError):
        rs.decode(surviving, len(data))


def test_decode_works_with_any_surviving_subset():
    """Erasure coding must recover regardless of *which* k blocks survive,
    not just a lucky subset."""
    data = os.urandom(500)
    rs = ReedSolomonErasureCoder(data_shards=8, parity_shards=4)
    blocks = rs.encode(data)
    for _ in range(10):
        random.shuffle(blocks)
        surviving = blocks[:8]
        decoded = rs.decode(surviving, len(data))
        assert decoded == data


def test_invalid_construction():
    with pytest.raises(ValueError):
        ReedSolomonErasureCoder(data_shards=0, parity_shards=2)
    with pytest.raises(ValueError):
        ReedSolomonErasureCoder(data_shards=5, parity_shards=-1)
