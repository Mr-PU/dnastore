# Contributing to dnastore

Thanks for considering a contribution. This project is a research/demo
storage simulator, not a production system -- contributions that make
it more correct, better documented, or more useful for experimentation
are all welcome.

## Getting started

```bash
git clone <your-repo-url>
cd dnastore
pip install -e ".[dev,vfs]"
pytest -v
python smoke_test.py   # dependency-free fallback if pytest isn't available
```

## Where things are

See the "Project layout" section in README.md. Briefly:
- `dnastore/codec/` -- binary <-> DNA base encodings
- `dnastore/ecc/` -- Reed-Solomon erasure coding (GF(256), pure Python)
- `dnastore/synthesis/`, `dnastore/sequencing/` -- simulators + real-vendor stubs
- `dnastore/api.py` -- the public `DNAStorage` interface

## Good first contributions

- **Real vendor backends**: `synthesis/twist_api.py` and
  `sequencing/nanopore_api.py` are deliberately unimplemented stubs. If
  you have API/lab access, filling these in against a real vendor's
  actual API is probably the single most valuable contribution.
- **Performance**: the Reed-Solomon path is O(k^3) pure-Python matrix
  inversion (see README's Performance section) -- a numpy-vectorized
  GF(256) implementation would help a lot for larger files.
- **New codecs**: implement `dnastore.codec.base.Codec` and register it
  in `dnastore.codec.REGISTRY`. `benchmarks/codec_comparison.py` and
  `benchmarks/efficiency_report.py` will pick it up automatically for
  comparison once you add it to their codec dicts.
- **Better error models**: the synthesis/sequencing simulators use
  fairly simple independent per-base error rates; real error profiles
  are correlated (e.g. GC-content-dependent dropout, position-dependent
  substitution rates) and more realistic models would make the
  simulator more useful for research.

## Before submitting a PR

- Run `pytest -v` (or `python smoke_test.py` if you don't have pytest)
  and make sure everything passes.
- If you touch the fountain codec or Reed-Solomon coder, run with
  multiple random seeds / trials -- both have probabilistic behavior
  that a single lucky run can hide bugs in (this bit us during initial
  development -- see CHANGELOG for the MDS-property bug that only
  showed up at larger shard counts).
- Update CHANGELOG.md with a one-line summary.
- If your change affects performance, re-run and update the numbers in
  README's "Performance characteristics" table rather than leaving
  stale numbers in place.

## Reporting bugs

Please include: Python version, codec used, file size, and (if
relevant) the specific error-rate parameters you passed to
`synthesis_kwargs`/`sequencing_kwargs` -- both the fountain decoder and
the Reed-Solomon reconstruction path have made assumptions that turned
out to be wrong at certain scales before, so exact repro parameters
matter more than usual here.
