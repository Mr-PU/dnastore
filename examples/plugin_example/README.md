# dnastore-xor-demo -- example codec plugin

A complete, working example of a third-party dnastore codec plugin. The
codec itself (`XorDemoCodec`) is deliberately trivial -- it XORs bytes
with a fixed key before delegating to dnastore's built-in `NaiveCodec`
for the actual DNA mapping. The point isn't the algorithm; it's proving
the plugin mechanism works end to end with a real, separately
installed package.

## Try it

```bash
pip install -e .          # from this directory
python3 -c "
import dnastore
print(dnastore.available_codecs())   # now includes 'xor_demo'

from dnastore import DNAStorage
dna = DNAStorage(codec='xor_demo')
dna.store(b'hello', name='f.txt')
print(dna.retrieve('f.txt'))
"
```

No changes to the dnastore repository were needed for `xor_demo` to
show up in `available_codecs()` and work through the full `DNAStorage`
API -- that's the whole point of entry-point-based discovery.

## How it works

Two things make this a plugin:

1. **A class implementing `dnastore.codec.base.Codec`** (see `dnastore_xor_demo/codec.py`) -- just `encode()` and `decode()`.
2. **An entry point declaration** in `pyproject.toml`:

```toml
[project.entry-points."dnastore.codecs"]
xor_demo = "dnastore_xor_demo.codec:XorDemoCodec"
```

That's the entire contract. `dnastore.plugins.discover_plugin_codecs()` scans installed packages' metadata for the `dnastore.codecs` entry-point group via `importlib.metadata` -- the same mechanism pytest, mkdocs, and much of the scientific Python ecosystem use for plugins.

## Writing a real plugin

Copy this directory as a starting point. Things to actually think about for a non-toy codec:

- **`self_redundant`**: set `True` on your `Codec` subclass if your scheme provides its own cross-strand redundancy (like the built-in `fountain` codec). Leave `False` (the default) if DNAStorage should wrap your codec's output in Reed-Solomon erasure coding across strands -- true for most "one shard in, one strand out" schemes.
- **`_shard_capacity_bytes`**: if `self_redundant` is `False`, declare this class attribute (the max input bytes your codec reliably packs into exactly one DNA strand) so DNAStorage can size Reed-Solomon shards correctly. If you don't declare it, `shard_capacity_bytes()` auto-probes by calling your `encode()` with increasing input sizes -- works, but slower and only as correct as your codec's behavior at the probed sizes.
- **Package naming convention**: `dnastore-<something>` is the suggested convention (mirrors `pytest-<plugin>`), though nothing enforces it.
