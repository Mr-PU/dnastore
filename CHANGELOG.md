# Changelog

## 0.1.0 -- initial release

- FASTA/FASTQ export (`dnastore.formats`, plus `DNAStorage.export_object_fasta()`,
  `export_pool_fasta()`, `retrieve_as_fastq()`): connects the simulated
  strand pool to the actual bioinformatics ecosystem -- FASTA is what
  real synthesis vendors accept for order submission, FASTQ is what
  sequencers/basecallers emit. Includes optional Biopython `SeqRecord`
  interop (`pip install dnastore[bio]`).
- Codec layer: `naive` (2-bit baseline), `rotating` (Goldman-style ternary,
  homopolymer-free by construction), `fountain` (LT/robust-soliton rateless
  code with self-verifying redundancy).
- Reed-Solomon erasure coding across strands (GF(256), pure Python, no
  external dependencies) for `naive`/`rotating`, using a full Vandermonde
  generator matrix (not a systematic identity+Vandermonde construction) so
  the MDS property -- any k of n surviving shards are always invertible --
  actually holds for every erasure pattern, not just the ones tested by
  chance.
- Automatic Reed-Solomon striping for large files: GF(256) caps a single
  codeword at 255 total shards, so files needing more data shards than
  that are split into independent RS-coded stripes transparently.
- Primer-pair addressing for selective/random-access reads.
- Versioned, tombstone-based metadata store modeling DNA's append-only
  nature (no in-place update).
- `SynthesisSimulator` / `SequencingSimulator` modeling substitution,
  indel, and (dominant) strand-dropout error modes without lab hardware.
- Integration stubs for real synthesis (`TwistSynthesisBackend`) and
  sequencing (`NanoporeSequencingBackend`) vendors.
- Optional FUSE-based virtual filesystem mount (`dnastore.vfs`).
- `DNAStorage` public API: `store()`, `retrieve()`, `update()`, `delete()`,
  `purge()`, `list_versions()`, `list_objects()`, `stats()`.
- Test suite (pytest) + a pytest-free `smoke_test.py` for offline/no-network
  environments.
- Benchmark script comparing codecs on density, GC-content stability,
  homopolymer risk, and dropout tolerance.
- Dockerfile + docker-compose for tests, benchmarks, examples, and the
  VFS mount.
