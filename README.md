# Begum: Metabarcoding preprocessing tool

Begum demultiplexes pooled metabarcoding / eDNA FASTQ reads by primer and
tag sequences (**sort**), then filters the resulting amplicons across PCR
replicates (**filter**).  It is available in two implementations:

| | Python 3 (`begum/`) | Rust (`begum-rs/`) |
|---|---|---|
| Requirements | Python ≥ 3.8, Biopython ≥ 1.78, regex, textdistance | Rust stable (cargo) |
| Input | plain or gzip FASTQ | plain or gzip FASTQ |
| Paired-end | yes | yes |
| Merged paired-end | yes (treated as single-end) | yes (treated as single-end) |

**Merged paired-end reads** (produced by AdapterRemoval, FLASH, PANDA, etc.) are
passed as single-end input — one FASTQ per pool, two columns in `pools.txt`.
Begum distinguishes single-end from paired-end solely by the number of columns
in `pools.txt`.

## Performance

Benchmarked on 100,000 single-end reads (2 pools, 2 samples, 2 replicates):

| Step | Python 3 | Rust | Speedup |
|------|----------|------|---------|
| sort | 11.8 s | 0.29 s | ~40× |
| filter | 0.48 s | 0.004 s | ~120× |
| **total** | **12.3 s** | **0.30 s** | **~40×** |

## Quick start

### Python version

```bash
pip install biopython regex textdistance

python3 begum/Begum.py sort \
  -p primers.txt -t tags.txt -s samples.txt -l pools.txt \
  -d out/ -o run1

python3 begum/Begum.py filter \
  -i out/run1 -s samples.txt \
  -p 1.0 -m 2 -l 50 \
  -d out/ -o filtered
```

### Rust version

```bash
cd begum-rs && cargo build --release
# binary: begum-rs/target/release/begum

begum sort \
  -p primers.txt -t tags.txt -s samples.txt -l pools.txt \
  -d out/ -o run1

begum filter \
  -i out/run1 -s samples.txt \
  -p 1.0 -m 2 -l 50 \
  -d out/ -o filtered
```

## Development history

Begum was originally written in Python 2.  [Claude Code](https://claude.ai/code)
was used to modernise and extend the codebase across several sessions:

1. **Python 2 → Python 3 port.**  Three compatibility issues were fixed:
   the `version=` keyword removed from `argparse.ArgumentParser` in Python 3;
   `Bio.Alphabet` and its IUPAC classes removed in Biopython ≥ 1.78; and
   `gzip.open()` returning bytes rather than text in Python 3.

2. **End-to-end tests.**  A Python test suite (`test/run_test.py`) was written
   covering single-end and paired-end sort and filter, with synthetic gzipped
   FASTQ inputs generated at runtime.

3. **Rust rewrite.**  Claude Code wrote an implementation plan
   (`docs/plans/2026-03-06-rust-port.md`) then executed it using the
   [Superpowers plugin](https://github.com/superpowers-sh/superpowers) for
   Claude Code, specifically the `subagent-driven-development` skill — a
   fresh subagent per task with spec-compliance and code-quality review after
   each.  The result is
   `begum-rs/`: a full reimplementation using
   [needletail](https://github.com/onecodex/needletail) for FASTQ parsing,
   [clap](https://github.com/clap-rs/clap) for the CLI, and a sliding-window
   IUPAC primer matcher in place of the Python `regex` fuzzy-match approach.
   The Rust version is ~40× faster than the Python version on a 100,000-read
   dataset.

4. **Tutorial and dataset.**  A 100,000-read synthetic dataset
   (`test/generate_complex_dataset.py`) was generated to demonstrate every
   sort outcome type (C/B/F/R/N, IUPAC primer matching, reverse-orientation
   reads) and every filter outcome (pass, fail propPCRs, fail minOccurrence,
   fail minLength).

5. **Repository cleanup.**  The repo was restructured to its current
   three-directory layout (`begum/`, `begum-rs/`, `test/`), a root
   `.gitignore` was added, and a pre-existing bug in the Python `filter`
   default for `--minLength` was fixed.

## Documentation

See `test/tutorial/README.md` for a full walkthrough covering all input file
formats, every command-line argument, output file formats, and a 100,000-read
example dataset that demonstrates all sort and filter outcomes.

## Testing

```bash
# Python end-to-end test
python3 test/run_test.py

# Rust unit + integration tests
cargo test --manifest-path begum-rs/Cargo.toml
```

## Repository layout

```
begum/       Python 3 source
begum-rs/    Rust source
test/        Tests and tutorial dataset
```
