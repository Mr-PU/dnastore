import os
import pytest

from dnastore import DNAStorage, DNAStorageError


@pytest.fixture(params=["naive", "rotating", "fountain"])
def dna(request):
    return DNAStorage(codec=request.param)


def test_store_and_retrieve_bytes(dna):
    data = os.urandom(2000)
    dna.store(data, name="file.bin")
    out = dna.retrieve("file.bin")
    assert out == data


def test_store_and_retrieve_file(dna, tmp_path):
    path = tmp_path / "doc.txt"
    path.write_bytes(b"hello dna storage" * 50)
    dna.store(str(path))
    out = dna.retrieve("doc.txt")
    assert out == path.read_bytes()


def test_update_creates_new_version(dna):
    dna.store(b"version 1", name="f.txt")
    dna.update(b"version 2", name="f.txt")
    versions = dna.list_versions("f.txt")
    assert len(versions) == 2
    assert dna.retrieve("f.txt") == b"version 2"


def test_delete_is_tombstone_not_erase(dna):
    dna.store(b"data", name="f.txt")
    dna.delete("f.txt")
    assert "f.txt" not in dna.list_objects()
    assert "f.txt" in dna.list_objects(include_tombstoned=True)
    with pytest.raises(FileNotFoundError):
        dna.retrieve("f.txt")


def test_purge_removes_physical_strands(dna):
    dna.store(b"data", name="f.txt")
    stats_before = dna.stats()
    assert stats_before["total_strands"] > 0
    dna.purge("f.txt")
    stats_after = dna.stats()
    assert stats_after["total_strands"] < stats_before["total_strands"]
    assert "f.txt" not in dna.list_objects(include_tombstoned=True)


def test_retrieve_missing_raises(dna):
    with pytest.raises(FileNotFoundError):
        dna.retrieve("does-not-exist.txt")


def test_state_persistence(tmp_path):
    state_path = str(tmp_path / "state.json")
    pool_path = str(tmp_path / "pool.json")

    dna1 = DNAStorage(codec="rotating", state_path=state_path, pool_path=pool_path)
    dna1.store(b"persisted content", name="f.txt")

    dna2 = DNAStorage(codec="rotating", state_path=state_path, pool_path=pool_path)
    assert dna2.retrieve("f.txt") == b"persisted content"


def test_store_requires_name_for_bytes():
    dna = DNAStorage()
    with pytest.raises(ValueError):
        dna.store(b"data")
