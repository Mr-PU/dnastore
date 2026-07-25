"""
Metadata store for the DNA object archive.

DNA cannot be rewritten in place, so this behaves like append-only
object storage with versioning (think S3 object versions or an
LSM-tree), not like a traditional in-place-updatable filesystem:

- store() creates version 1 of an object.
- update() synthesizes a brand new strand set and adds a new version;
  the old version's strands are *not* physically erased, only marked
  obsolete in metadata.
- delete() writes a tombstone marker; it does not free any DNA.
- purge() is the explicit, rarer operation that actually discards a
  version's underlying strand records -- modeling a real physical
  "destroy this sample" operation.
"""
from __future__ import annotations

import time
import uuid
from dataclasses import dataclass, field


@dataclass
class Version:
    version_id: str
    version_number: int
    created_at: float
    size_bytes: int
    checksum: str
    codec: str
    ecc: str | None
    strand_count: int
    obsolete: bool = False


@dataclass
class ObjectRecord:
    object_id: str
    name: str
    created_at: float
    versions: list[Version] = field(default_factory=list)
    tombstoned: bool = False
    tombstoned_at: float | None = None

    def latest_version(self) -> Version | None:
        active = [v for v in self.versions if not v.obsolete]
        return active[-1] if active else None


class MetadataStore:
    def __init__(self):
        self._objects: dict[str, ObjectRecord] = {}
        self._name_to_id: dict[str, str] = {}

    def create_object(self, name: str) -> ObjectRecord:
        if name in self._name_to_id:
            existing = self._objects[self._name_to_id[name]]
            if not existing.tombstoned:
                raise FileExistsError(
                    f"object '{name}' already exists; use update() to create a new version"
                )
        object_id = str(uuid.uuid4())
        record = ObjectRecord(object_id=object_id, name=name, created_at=time.time())
        self._objects[object_id] = record
        self._name_to_id[name] = object_id
        return record

    def add_version(
        self, name: str, size_bytes: int, checksum: str, codec: str,
        ecc: str | None, strand_count: int, mark_previous_obsolete: bool = True,
    ) -> Version:
        record = self.get_by_name(name, include_tombstoned=True)
        if record is None:
            record = self.create_object(name)
        else:
            record.tombstoned = False
            record.tombstoned_at = None

        if mark_previous_obsolete:
            for v in record.versions:
                v.obsolete = True

        version = Version(
            version_id=str(uuid.uuid4()),
            version_number=len(record.versions) + 1,
            created_at=time.time(),
            size_bytes=size_bytes,
            checksum=checksum,
            codec=codec,
            ecc=ecc,
            strand_count=strand_count,
            obsolete=False,
        )
        record.versions.append(version)
        return version

    def get_by_name(self, name: str, include_tombstoned: bool = False) -> ObjectRecord | None:
        object_id = self._name_to_id.get(name)
        if object_id is None:
            return None
        record = self._objects[object_id]
        if record.tombstoned and not include_tombstoned:
            return None
        return record

    def tombstone(self, name: str) -> None:
        record = self.get_by_name(name)
        if record is None:
            raise KeyError(f"no such object: '{name}'")
        record.tombstoned = True
        record.tombstoned_at = time.time()

    def purge(self, name: str) -> ObjectRecord:
        """Remove the object's metadata entirely. Returns the removed
        record so the caller can also purge underlying strand records
        from the physical pool."""
        object_id = self._name_to_id.pop(name, None)
        if object_id is None:
            raise KeyError(f"no such object: '{name}'")
        return self._objects.pop(object_id)

    def list_objects(self, include_tombstoned: bool = False) -> list[ObjectRecord]:
        return [
            r for r in self._objects.values()
            if include_tombstoned or not r.tombstoned
        ]

    def to_dict(self) -> dict:
        def version_dict(v: Version) -> dict:
            return {
                "version_id": v.version_id, "version_number": v.version_number,
                "created_at": v.created_at, "size_bytes": v.size_bytes,
                "checksum": v.checksum, "codec": v.codec, "ecc": v.ecc,
                "strand_count": v.strand_count, "obsolete": v.obsolete,
            }

        return {
            "objects": {
                oid: {
                    "object_id": r.object_id, "name": r.name, "created_at": r.created_at,
                    "tombstoned": r.tombstoned, "tombstoned_at": r.tombstoned_at,
                    "versions": [version_dict(v) for v in r.versions],
                }
                for oid, r in self._objects.items()
            },
            "name_to_id": dict(self._name_to_id),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "MetadataStore":
        store = cls()
        for oid, obj in data.get("objects", {}).items():
            record = ObjectRecord(
                object_id=obj["object_id"], name=obj["name"], created_at=obj["created_at"],
                tombstoned=obj.get("tombstoned", False), tombstoned_at=obj.get("tombstoned_at"),
            )
            for v in obj.get("versions", []):
                record.versions.append(Version(**v))
            store._objects[oid] = record
        store._name_to_id = dict(data.get("name_to_id", {}))
        return store
