# dnastore

[![CI](https://github.com/yourname/dnastore/actions/workflows/ci.yml/badge.svg)](https://github.com/yourname/dnastore/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/dnastore.svg)](https://pypi.org/project/dnastore/)

A programmable storage abstraction over DNA. `store()`, `retrieve()`, `update()`, `delete()` -- same shape as an object storage SDK, but backed by a simulated (or, with real synthesis/sequencing hardware, actual) DNA archive.

```python
from dnastore import DNAStorage

dna = DNAStorage()
dna.store("report.pdf")
data = dna.retrieve("report.pdf")
dna.update("report.pdf")   # new version; old strands marked obsolete, not erased
dna.delete("report.pdf")   # tombstone; DNA is not physically freed until purge()
```

## Why this exists

DNA can't be rewritten in place. Once a strand is synthesized, "updating" it means synthesizing a new strand and marking the old one obsolete -- which makes DNA storage behave much more like **append-only object storage** (think S3 versioning, or an LSM-tree) than like a disk. This project takes that seriously as the core design constraint, rather than pretending DNA is a block device.

It's runnable **today**, with zero lab access, via a physically-grounded simulator that injects the error modes real synthesis and sequencing actually produce (substitution errors, indels, and -- the dominant failure mode -- whole-strand dropout). The same API has a clean seam for real synthesis/sequencing vendor integration later.

## Architecture

```
Application
      |
      v
DNAStorage API           store() / retrieve() / delete() / update()
      |
      v
Codec                     binary <-> DNA bases (naive / rotating / fountain)
      |
      v
Error Correction           Reed-Solomon erasure coding across strands
      |                    (fountain codec is self-redundant, skips this layer)
      v
Addressing                 primer-pair allocation (selective/random access)
      |
      v
Synthesis                  writes strands (simulator, or a real vendor backend)
      |
      v
Physical Pool               the "physical DNA" -- in-memory/on-disk archive
      |
      v
Sequencing                  reads strands back (simulator, or a real backend)
      |
      v
Codec decode  -->  Original Data
```

## Codecs

| Codec | Bits/base | Homopolymer risk | Redundancy strategy |
|---|---|---|---|
| `naive` | 2.0 (max density) | High -- long runs possible | External Reed-Solomon erasure coding |
| `rotating` | ~1.58 | None by construction (Goldman-style ternary encoding) | External Reed-Solomon erasure coding |
| `fountain` | Adaptive (~0.5-1.5 depending on loss tolerance) | Same as naive per-droplet | Self-redundant (rateless LT/fountain code, robust soliton distribution) |

Run `python benchmarks/codec_comparison.py` for a live comparison table (density, GC-content stability, homopolymer runs, empirical dropout tolerance) on your machine.

`fountain` is the default: it's what real DNA storage systems (e.g. the "DNA Fountain" architecture from Erlich & Zielinski, 2017) use for exactly this reason -- rateless codes tolerate strand dropout gracefully instead of needing a separately-tuned erasure code.

## Plugin architecture

Third-party packages can register new codecs without touching this repository, via Python entry points -- the same mechanism pytest and mkdocs use for plugins:

```python
import dnastore
print(dnastore.available_codecs())   # ['fountain', 'naive', 'rotating', ...plus any installed plugins]

dna = dnastore.DNAStorage(codec="some_plugin_codec")
```

A plugin package declares itself in its own `pyproject.toml`:

```toml
[project.entry-points."dnastore.codecs"]
my_codec = "my_package.module:MyCodecClass"
```

Once `pip install`-ed alongside dnastore, `my_codec` shows up in `available_codecs()` and works through the full `DNAStorage` API automatically -- no PR against this repo required. See `examples/plugin_example/` for a complete, real working example (a full separate installable package with its own `pyproject.toml`, actually built and installed as part of this project's own test suite to prove the mechanism works, not just described).

## Why DNA storage needs its own error model

Unlike an SSD, DNA's dominant failure mode isn't a flipped bit -- it's a strand that never comes back at all (failed synthesis, degraded storage, failed PCR amplification). This project models that explicitly:

- Every strand gets a CRC16 checksum. A strand that fails its checksum on readback is treated as **erased**, not corrected in place -- correction is delegated entirely to the cross-strand layer (Reed-Solomon or the fountain decoder's own redundancy).
- `SynthesisSimulator` injects per-base substitution errors and per-strand dropout.
- `SequencingSimulator` injects substitution errors, insertions/deletions, and read-time dropout.
- Default error rates are tuned to be realistic-but-comfortable (~10-15% expected per-strand loss after both stages); crank them up via `synthesis_kwargs`/`sequencing_kwargs` to see where your redundancy budget breaks down -- see `examples/stress_test.py`.

## Bioinformatics interop (FASTA/FASTQ)

FASTA and FASTQ are the formats that actually connect this project to the rest of the bioinformatics world: every real synthesis vendor accepts FASTA order submissions, and every sequencer/basecalling pipeline emits FASTQ. `dnastore.formats` implements both natively (no dependencies), plus optional Biopython `SeqRecord` conversion:

```python
dna.export_object_fasta("report.pdf", "synthesis_order.fasta")  # what you'd submit to a vendor
dna.export_pool_fasta("entire_archive.fasta")                    # whole physical pool, all objects
dna.retrieve_as_fastq("report.pdf", "reads.fastq")                # raw sequencer-style reads + quality scores

# optional: pip install dnastore[bio]
from dnastore.formats import to_seqrecord, read_fasta
records = [to_seqrecord(r) for r in read_fasta("synthesis_order.fasta")]  # -> Bio.SeqRecord objects
```

See `examples/bio_interop.py` for a full walkthrough.

## Random access via primer addressing

Real DNA storage doesn't sequence an entire archive to read one file -- it flanks each object's strands with a unique PCR primer pair and selectively amplifies just those strands. `dnastore.addressing.PrimerRegistry` models this: every object gets a unique primer pair, and `SequencingBackend.sequence(primer_fwd, primer_rev)` only returns strands for that object. This is what makes `retrieve()` a targeted read instead of "decode everything in the pool."

## Versioning and delete semantics

- `store()` / `update()` create a new version; the previous version's strands are marked obsolete in metadata but not physically removed (DNA is append-only).
- `delete()` writes a tombstone -- the object disappears from `list_objects()` and `retrieve()` raises `FileNotFoundError`, but nothing is physically erased.
- `purge()` is the explicit, rarer, irreversible operation that actually discards the underlying strand records -- modeling a real "destroy this sample" lab operation.

## Installation

```bash
pip install -e .                 # core library
pip install -e ".[vfs]"          # + FUSE mount support
pip install -e ".[bio]"          # + Biopython SeqRecord interop
pip install -e ".[dev]"          # + pytest
```

Or via Docker (see below) -- no local Python setup needed.

## Quickstart

```bash
python examples/quickstart.py
python examples/stress_test.py       # error-rate tuning demo
python benchmarks/codec_comparison.py
```

```python
from dnastore import DNAStorage

dna = DNAStorage(
    codec="fountain",                        # or "naive", "rotating"
    state_path="archive_state.json",         # persists metadata/versions across restarts
    pool_path="archive_pool.json",           # persists the physical strand pool
)

dna.store(b"hello world", name="hello.txt")  # or dna.store("path/to/file.pdf")
data = dna.retrieve("hello.txt")

dna.update(b"new content", name="hello.txt")
print(dna.list_versions("hello.txt"))

dna.delete("hello.txt")     # tombstone
dna.purge("hello.txt")      # actually free the DNA
```

## Mounting as a real filesystem (FUSE)

```bash
pip install -e ".[vfs]"          # requires FUSE installed at the OS level too
python examples/mount_vfs.py /mnt/dna

# in another terminal:
cp somefile.txt /mnt/dna/
cat /mnt/dna/somefile.txt
rm /mnt/dna/somefile.txt
```

Writes are buffered per file descriptor and only actually synthesized on `close()`, since DNA can't be appended to incrementally.

## Docker

```bash
docker build -t dnastore .
docker run --rm dnastore                          # smoke test
docker run --rm dnastore python examples/quickstart.py
docker run --rm dnastore python benchmarks/codec_comparison.py

# or via docker-compose:
docker compose run --rm test
docker compose run --rm benchmark
docker compose run --rm quickstart

# FUSE mount needs elevated container permissions:
docker compose run --rm vfs
```

## Real hardware integration (future)

`dnastore/synthesis/twist_api.py` and `dnastore/sequencing/nanopore_api.py` are integration stubs -- not working clients. Real synthesis is an async, multi-day-turnaround process gated by vendor account access, and real sequencing needs a basecalling/demultiplexing pipeline in front of it. They exist to mark where that integration goes:

```python
dna = DNAStorage(backend="twist+nanopore", synthesis_kwargs={"api_key": "..."})
```

If you have lab or vendor API access and want to build this out, contributions welcome -- see the stub files for exactly what's missing.

## Project layout

```
dnastore/
  api.py              # DNAStorage -- the public interface
  codec/               # naive / rotating / fountain codecs
  ecc/                 # Reed-Solomon erasure coding (GF(256))
  addressing.py         # primer-pair registry
  metadata.py           # versioning, tombstones
  store.py              # PhysicalPool -- the "physical DNA" abstraction
  synthesis/            # SynthesisSimulator + real-vendor stub
  sequencing/           # SequencingSimulator + real-vendor stub
  vfs.py                # optional FUSE mount
  checksum.py           # CRC16 for per-strand erasure detection
  formats.py             # FASTA/FASTQ export + optional Biopython interop
  gf256.py              # Galois field arithmetic backing Reed-Solomon
tests/                  # pytest suite
smoke_test.py           # dependency-free end-to-end check (no pytest needed)
benchmarks/              # codec comparison
examples/                # quickstart, VFS mount, error-rate stress test
```

## Testing

```bash
pytest                    # full suite (needs pip install -e ".[dev]")
python smoke_test.py       # no dependencies beyond the stdlib + this package
```

## Is DNA storage actually "more efficient"?

Not in the way that usually matters day-to-day. Speed and cost are DNA storage's weakest points today, by a wide margin -- this project's own numbers make that obvious (see Performance below, and `benchmarks/efficiency_report.py`). Where DNA genuinely wins is **physical density and durability**: roughly 10 orders of magnitude denser than a hard drive by mass, and stable for centuries under the right conditions, versus decades for magnetic media.

`benchmarks/efficiency_report.py` measures each codec's *actual* bases-per-byte through the real `DNAStorage` pipeline (including CRC and redundancy overhead), then converts that into real mass (from DNA's molecular weight -- ~330 g/mol per nucleotide, first-principles physics, not a borrowed headline number) and cost (using published per-base synthesis pricing). Run it yourself:

```bash
python benchmarks/efficiency_report.py
```

The result: codec choice changes cost/mass by roughly 1.5-3x (redundancy overhead). Synthesis cost versus an SSD today differs by a factor of 10^5-10^8. The codec is a rounding error next to the actual bottleneck.

## Performance characteristics (measured, not estimated)

This is a pure-Python implementation with no numpy/C acceleration in the
hot paths, so scaling behavior differs meaningfully by codec:

| File size | `fountain` (store/retrieve) | `rotating`+RS (store/retrieve) | `naive`+RS (store/retrieve) |
|---|---|---|---|
| 50 KB | 0.7s / 0.7s | 1.4s / 6.3s | -- |
| 200 KB | 2.8s / 3.1s | 5.5s / 26s | -- |
| 500 KB | 7.1s / 8.2s | -- | 12s / 45s |

`fountain` scales roughly linearly (XOR-based peeling decode). The
Reed-Solomon-backed codecs (`naive`, `rotating`) scale worse because
GF(256) matrix inversion during decode is O(k³) per stripe in pure
Python -- fine for small files and for understanding the algorithm,
but genuinely slow for anything beyond a few hundred KB. **For files
above a few hundred KB, prefer `codec="fountain"`**, or treat the
RS path as a reference implementation to optimize (numpy-vectorized
GF(256) arithmetic, or swapping in a C-accelerated erasure-coding
library) rather than something to run as-is at scale.

## Known limitations / honest caveats

- The Reed-Solomon layer corrects **erasures** (whole missing/corrupted strands, detected via CRC), not in-strand substitution errors directly -- that's a deliberate design choice matching how real systems separate "did this strand survive intact" from "reconstruct what's missing," but it means a single bad base anywhere in a strand costs you the whole strand rather than just that base.
- GF(256) caps a single Reed-Solomon codeword at 255 total shards. Files needing more data shards than that are automatically split into independent RS-coded stripes (`DNAStorage.STRIPE_MAX_DATA_SHARDS`, default 150 data shards/stripe) -- the same principle as RAID striping across multiple parity groups. This is handled transparently; you don't need to think about it unless you're tuning very large archives.
- The fountain codec's decoder is a from-scratch implementation of the Luby Transform with a robust soliton distribution; it self-verifies decodability at encode time (growing redundancy until an internal check passes) rather than relying purely on closed-form overhead formulas, which is a pragmatic engineering choice, not a claim of asymptotically optimal overhead.
- `MAX_STRAND_LENGTH` (200 nt) and shard-capacity constants are illustrative defaults based on typical synthesis vendor limits, not tied to any specific vendor's actual current specs -- confirm against your vendor before relying on them.
- This is a research/demo-grade implementation. It has not been audited for production use, and the real-hardware backends are unimplemented stubs.

## Contributing

See `CONTRIBUTING.md`.

## Releasing (for maintainers)

CI (`.github/workflows/ci.yml`) runs on every push/PR: pytest across Python 3.10-3.12, the smoke test, and a Docker build+run.

To publish a new version to PyPI:

1. Bump the version in `pyproject.toml` and `dnastore/__init__.py` (`__version__`), update `CHANGELOG.md`.
2. Configure [PyPI trusted publishing](https://docs.pypi.org/trusted-publishers/) for this repo once (Settings on the PyPI project page -- no long-lived API token needed after that).
3. `git tag v0.1.1 && git push origin v0.1.1` -- this triggers `.github/workflows/publish.yml`, which builds, validates, smoke-tests the built wheel, and publishes automatically.

Or manually, without the workflow:
```bash
pip install build twine
python -m build
twine check dist/*
twine upload --repository testpypi dist/*   # dry run first
twine upload dist/*                          # then the real thing
```

## License

Apache 2.0 -- see `LICENSE`. If you build on the DNA Fountain (Erlich & Zielinski, 2017) or Goldman et al. (2013) rotating-code ideas, please cite the original papers.
