"""
Virtual filesystem mount for a DNAStorage archive, using FUSE.

    dna = DNAStorage(state_path="archive_state.json", pool_path="archive_pool.json")
    mount_dna_fs(dna, "/mnt/dna")

Then ordinary tools work directly against DNA-backed storage:

    cp report.pdf /mnt/dna/
    cat /mnt/dna/report.pdf > restored.pdf
    rm /mnt/dna/report.pdf   # tombstones, doesn't physically erase

This is a read-mostly, whole-file filesystem: writes are buffered in
memory per file descriptor and only committed to DNAStorage (i.e.
actually synthesized) on close(), since DNA can't be updated
in place or written incrementally the way a real block device can.

Requires the optional `fusepy` dependency and a working FUSE
installation on the host (`apt-get install fuse` on Debian/Ubuntu, or
macFUSE on macOS). Install with: pip install dnastore[vfs]
"""
from __future__ import annotations

import errno
import stat
import time

try:
    import fuse
    from fuse import FUSE, FuseOSError, Operations
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "The VFS mount requires the optional 'fusepy' dependency and a "
        "working FUSE installation. Install with: pip install dnastore[vfs] "
        "(and on Linux, ensure the 'fuse' package is installed at the OS level)."
    ) from e

from .api import DNAStorage


class DNAFS(Operations):
    def __init__(self, dna: DNAStorage):
        self.dna = dna
        self._write_buffers: dict[str, bytearray] = {}
        self._mount_time = time.time()

    # -- helpers --------------------------------------------------------

    def _strip(self, path: str) -> str:
        return path.lstrip("/")

    def _exists(self, name: str) -> bool:
        return name in self.dna.list_objects()

    # -- filesystem metadata ---------------------------------------------

    def getattr(self, path, fh=None):
        if path == "/":
            return {
                "st_mode": stat.S_IFDIR | 0o755, "st_nlink": 2,
                "st_ctime": self._mount_time, "st_mtime": self._mount_time, "st_atime": self._mount_time,
            }
        name = self._strip(path)
        if name in self._write_buffers:
            size = len(self._write_buffers[name])
        elif self._exists(name):
            version = self.dna.metadata.get_by_name(name).latest_version()
            size = version.size_bytes
        else:
            raise FuseOSError(errno.ENOENT)
        return {
            "st_mode": stat.S_IFREG | 0o644, "st_nlink": 1, "st_size": size,
            "st_ctime": self._mount_time, "st_mtime": self._mount_time, "st_atime": self._mount_time,
        }

    def readdir(self, path, fh):
        return [".", ".."] + self.dna.list_objects()

    # -- file operations --------------------------------------------------

    def open(self, path, flags):
        return 0  # stateless; buffering keyed by name, not a real fd table

    def create(self, path, mode, fi=None):
        name = self._strip(path)
        self._write_buffers[name] = bytearray()
        return 0

    def read(self, path, size, offset, fh):
        name = self._strip(path)
        data = self._write_buffers.get(name)
        if data is None:
            data = self.dna.retrieve(name)
        return bytes(data[offset:offset + size])

    def write(self, path, data, offset, fh):
        name = self._strip(path)
        buf = self._write_buffers.setdefault(name, bytearray())
        if len(buf) < offset:
            buf.extend(b"\x00" * (offset - len(buf)))
        buf[offset:offset + len(data)] = data
        return len(data)

    def truncate(self, path, length, fh=None):
        name = self._strip(path)
        buf = self._write_buffers.setdefault(name, bytearray())
        del buf[length:]

    def flush(self, path, fh=None):
        return 0

    def release(self, path, fh):
        """This is where the actual DNA synthesis happens -- on close(),
        not incrementally, since DNA can't be appended to in place."""
        name = self._strip(path)
        if name in self._write_buffers:
            data = bytes(self._write_buffers.pop(name))
            if self._exists(name):
                self.dna.update(data, name=name)
            else:
                self.dna.store(data, name=name)
        return 0

    def unlink(self, path):
        name = self._strip(path)
        self.dna.delete(name)

    def rename(self, old, new):
        old_name, new_name = self._strip(old), self._strip(new)
        data = self.dna.retrieve(old_name)
        self.dna.store(data, name=new_name)
        self.dna.delete(old_name)


def mount_dna_fs(dna: DNAStorage, mountpoint: str, foreground: bool = True) -> None:
    """Blocks until the filesystem is unmounted (Ctrl-C or `fusermount -u`)."""
    FUSE(DNAFS(dna), mountpoint, foreground=foreground, nothreads=True)
