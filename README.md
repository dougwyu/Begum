# Begum: Metabarcoding preprocessing tool

Begum demultiplexes pooled metabarcoding / eDNA FASTQ reads by primer and
tag sequences (**sort**), then filters the resulting amplicons across PCR
replicates (**filter**).  It is available in two implementations:

| | Python 3 (`begum/`) | Rust (`begum-rs/`) |
|---|---|---|
| Requirements | Python ≥ 3.8, Biopython ≥ 1.78, regex, textdistance | Rust stable (cargo) |
| Input | plain or gzip FASTQ | plain or gzip FASTQ |
| Paired-end | yes | yes |

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
