"""Mount a DNAStorage archive as a real filesystem via FUSE.

Requires: pip install dnastore[vfs], plus a working FUSE install
(`apt-get install fuse` on Debian/Ubuntu, macFUSE on macOS).

Run: python examples/mount_vfs.py /path/to/mountpoint
Then, in another terminal:
    cp somefile.txt /path/to/mountpoint/
    cat /path/to/mountpoint/somefile.txt
    rm /path/to/mountpoint/somefile.txt
Unmount with Ctrl-C here, or `fusermount -u /path/to/mountpoint`.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dnastore import DNAStorage
from dnastore.vfs import mount_dna_fs


def main():
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <mountpoint>")
        sys.exit(1)

    mountpoint = sys.argv[1]
    os.makedirs(mountpoint, exist_ok=True)

    dna = DNAStorage(
        codec="rotating",
        state_path="/tmp/dnastore_vfs_state.json",
        pool_path="/tmp/dnastore_vfs_pool.json",
    )

    print(f"Mounting DNA-backed filesystem at {mountpoint} (Ctrl-C to unmount)...")
    mount_dna_fs(dna, mountpoint, foreground=True)


if __name__ == "__main__":
    main()
