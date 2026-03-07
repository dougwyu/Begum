# Begum Tutorial

Begum is a metabarcoding / eDNA preprocessing tool that demultiplexes FASTQ
reads by primer and tag sequences (the *sort* step), then filters the
resulting amplicons by quality criteria across PCR replicates (the *filter*
step).

This tutorial covers:
1. [Installing the Python version](#1-installation--python-version)
2. [Installing the Rust version](#2-installation--rust-version)
3. [Input file formats](#3-input-file-formats)
4. [Running `sort`](#4-running-sort)
5. [Sort output files](#5-sort-output-files)
6. [Running `filter`](#6-running-filter)
7. [Filter output file](#7-filter-output-file)
8. [Complex dataset walkthrough](#8-complex-dataset-walkthrough)

---

## 1. Installation — Python version

**Requirements:** Python 3.8+, Biopython ≥ 1.78, regex, textdistance.

```bash
# Using pip
pip install biopython regex textdistance

# Verify
python3 begum/Begum.py --version
# 0.1
```

All commands are run from the repository root:

```bash
python3 begum/Begum.py sort   [options]
python3 begum/Begum.py filter [options]
```

---

## 2. Installation — Rust version

**Requirements:** Rust toolchain (stable), `cargo`.  Install from
[rustup.rs](https://rustup.rs) if not present.

```bash
cd begum-rs
cargo build --release
# Binary is at: begum-rs/target/release/begum
```

Run from anywhere:

```bash
./begum-rs/target/release/begum sort   [options]
./begum-rs/target/release/begum filter [options]
```

Or install globally:

```bash
cargo install --path begum-rs
begum sort   [options]
begum filter [options]
```

---

## 3. Input file formats

### primers.txt

One line with two sequences separated by whitespace: forward primer and
reverse primer.  IUPAC ambiguity codes are supported (e.g. `R` = A or G).

```
GCRTGC AGTCAG
```

Alternatively, supply primers directly on the command line with `--p1` /
`--p2` instead of a file.

### tags.txt

One tag per line: tag name and its sequence (forward tag only — the reverse
tag is the reverse complement of the same tag sequence).

```
tag1 AACCGGT
tag2 TTGGCCA
tag3 CCGGAAT
tag4 GGCCTTA
tag5 AATCCGG
tag6 TTAAGGC
tag7 GGTAACC
tag8 CCAATTG
```

Tags should have a pairwise Hamming distance of at least 3 to avoid
ambiguous matching when mismatches are allowed.

### samples.txt

One PCR replicate per line: `SampleName FwdTagName RevTagName PoolName`.

```
Sample1 tag1 tag2 Pool1
Sample1 tag3 tag4 Pool2
Sample2 tag5 tag6 Pool1
Sample2 tag7 tag8 Pool2
```

- The same sample appearing more than once defines PCR replicates.  The
  first appearance is replicate 1, the second is replicate 2, etc.
- Tag combinations must be unique within each pool.

### pools.txt

One pool per line: `PoolName /path/to/reads.fastq.gz [/path/to/reads2.fastq.gz]`.

Single-end example (gzip-compressed):
```
Pool1 /data/reads_pool1.fastq.gz
Pool2 /data/reads_pool2.fastq.gz
```

Single-end example (uncompressed):
```
Pool1 /data/reads_pool1.fastq
Pool2 /data/reads_pool2.fastq
```

Paired-end example (both files required on every line):
```
Pool1 /data/pool1_R1.fastq.gz /data/pool1_R2.fastq.gz
```

- All pools must be either all single-end or all paired-end.
- Gzip compression is detected automatically from the `.gz` file extension.
- Begum detects single-end vs paired-end solely from the number of columns
  in `pools.txt` (2 = single-end, 3 = paired-end).

**Merged paired-end reads** (produced by AdapterRemoval, FLASH, PANDA, etc.)
are treated as single-end: supply one FASTQ per pool with two columns.  The
merged read has the same structure as a single-end read:

```
[fwd_tag][fwd_primer]───── amplicon ─────[rc(rev_primer)][rc(rev_tag)]
```

### FASTQ files

Standard FASTQ format.  For single-end (including merged paired-end) reads,
the expected structure is:

```
[fwd_tag][fwd_primer][amplicon][rc(rev_primer)][rc(rev_tag)]
```

Begum also handles reverse-orientation reads (e.g. some merging tools output
reads in either orientation):

```
[rev_tag][rev_primer][rc(amplicon)][rc(fwd_primer)][rc(fwd_tag)]
```

For unmerged paired-end reads the two files are sequenced from opposite ends
of the fragment, so each read only contains one primer and one tag:

```
Read 1:  [fwd_tag][fwd_primer][amplicon_R1 →
Read 2:  [rev_tag][rev_primer][amplicon_R2 →
```

---

## 4. Running `sort`

### Python

```bash
python3 begum/Begum.py sort \
  -p primers.txt \          # primer file
  -t tags.txt \             # tag file
  -s samples.txt \          # sample information file
  -l pools.txt \            # pool information file
  -d out/ \                 # output directory (default: .)
  -o prefix                 # output file prefix (default: empty)
```

### Rust

```bash
begum sort \
  -p primers.txt \
  -t tags.txt \
  -s samples.txt \
  -l pools.txt \
  -d out/ \
  -o prefix
```

### All sort arguments

| Flag | Long form | Default | Description |
|------|-----------|---------|-------------|
| `-p` | `--primers` | — | Primer file (two-column: fwd rev) |
| `--p1` | `--fwdPrimer` | — | Forward primer sequence (alternative to -p) |
| `--p2` | `--revPrimer` | — | Reverse primer sequence (alternative to -p) |
| `-t` | `--tags` | **required** | Tag file |
| `-s` | `--sampleInfo` | **required** | Sample information file |
| `-l` | `--pool` | **required** | Pool information file |
| `-m` | `--allowMultiplePrimers` | false | Allow more than one primer match per read |
| `--pm` | `--primerMismatches` | 0 | Allowed mismatches in primers |
| `--tm` | `--tagMismatches` | 0 | Allowed mismatches in tags |
| `-d` | `--output_directory` | `.` | Output directory |
| `-o` | `--output_prefix` | `""` | Prefix for output files |

Either `-p` or both `--p1` and `--p2` must be provided.

---

## 5. Sort output files

For each pool named `<Pool>`, sort writes two files:

- `<prefix>_<Pool>.tagInfo`
- `<prefix>_<Pool>.summaryCounts`

### `<prefix>_<Pool>.tagInfo`

Tab-separated.  One row per unique (fwd_tag, rev_tag, amplicon) combination
observed in the pool.

**Single-end columns:**

| Column | Description |
|--------|-------------|
| `FTag` | Forward tag name matched at the start of the read |
| `RTag` | Reverse tag name matched at the start of the reverse-complemented read |
| `Seq` | Amplicon sequence (between primers, tags stripped) |
| `Count` | Number of reads with this exact tag + amplicon combination |
| `Type` | Classification code (see below) |

**Paired-end adds `FSeq` and `RSeq` in place of `Seq`.**

#### Type codes

| Code | Meaning |
|------|---------|
| `C` | **Correct** — both tags are in the pool and form the expected pair |
| `B` | **Both** tags present in pool, but not the expected combination |
| `F` | **Forward** tag in pool; reverse tag not in pool |
| `R` | **Reverse** tag in pool; forward tag not in pool |
| `N` | **Neither** tag belongs to this pool |

Only `Type=C` rows are used downstream by the filter step.

Example Pool1 tagInfo from the tutorial dataset:

```
FTag    RTag    Seq                                                             Count   Type
tag1    tag2    AAATTTGGGCCC...AAATTTGGGCCC  (60 nt)                          24000   C
tag1    tag2    GCTT                                                             5000   C
tag1    tag2    GGGCCCAAATTT...GGGCCCAAATTT  (60 nt)                           5000   C
tag1    tag2    TTTGGGCCCAAA...TTTGGGCCCAAA  (60 nt)                               1   C
tag1    tag3    AAATTTGGGCCC...  (60 nt)                                        2000   F
tag5    tag6    AAATTTGGGCCC...  (60 nt)                                        5000   C
tag5    tag6    GCTT                                                             1000   C
tag5    tag6    TTTGGGCCCAAA...  (60 nt)                                            1   C
tag1    tag6    AAATTTGGGCCC...  (60 nt)                                        2000   B
tag3    tag2    AAATTTGGGCCC...  (60 nt)                                        2000   R
tag3    tag4    AAATTTGGGCCC...  (60 nt)                                        1000   N
```

Reads without any primer match do not appear in tagInfo at all.

### `<prefix>_<Pool>.summaryCounts`

Tab-separated summary of unique-sequence and total-read counts per tag
combination, grouped by type.

```
Tag1    Tag2    UniqSeqs    TotalSeqs   Type
Correct combination of tags used in pool
----------------------------------------
tag1    tag2    4           34001       C
tag5    tag6    3           6001        C
Both tags used in pool, but not this combination
------------------------------------------------
tag1    tag6    1           2000        B
Only forward tag used in pool, reverse not found
------------------------------------------------
tag1    tag3    1           2000        F
Only reverse tag used in pool, forward not found
------------------------------------------------
tag3    tag2    1           2000        R
Neither tag used in this pool
-----------------------------
tag3    tag4    1           1000        N
```

---

## 6. Running `filter`

### Python

```bash
python3 begum/Begum.py filter \
  -i out/prefix \          # prefix of .tagInfo files from sort (without _Pool*.tagInfo)
  -s samples.txt \         # same sample information file used in sort
  -p 1.0 \                 # propPCRs: minimum fraction of reps a sequence must appear in
  -m 2 \                   # minOccurrence: minimum read count per replicate
  -l 50 \                  # minLength: minimum amplicon length (nt)
  -d out/ \                # output directory
  -o filtered              # output file prefix
```

### Rust

```bash
begum filter \
  -i out/prefix \
  -s samples.txt \
  -p 1.0 \
  -m 2 \
  -l 50 \
  -d out/ \
  -o filtered
```

### All filter arguments

| Flag | Long form | Default | Description |
|------|-----------|---------|-------------|
| `-i` | `--inputPrefix` | **required** | Path prefix of sort output (before `_<Pool>.tagInfo`) |
| `-s` | `--sampleInfo` | **required** | Sample information file |
| `-p` | `--propPCRs` | `1.0` | Minimum proportion of replicates the sequence must pass in |
| `-m` | `--minOccurence` | `1` | Minimum read count in a replicate for it to count as "present" |
| `-l` | `--minLength` | `1` | Minimum amplicon length in nucleotides |
| `-d` | `--output_directory` | `.` | Output directory |
| `-o` | `--output_prefix` | `Filtered` | Output file prefix |

#### How propPCRs works

For each sample with *N* total replicates, a sequence passes the `propPCRs`
filter if:

```
(number of reps where count >= minOccurrence) / N  >=  propPCRs
```

With `propPCRs=1.0` the sequence must be present (above `minOccurrence`) in
**every** replicate.  With `propPCRs=0.5` it must appear in at least half.

---

## 7. Filter output file

### `<output_prefix>.fna`

FASTA format.  One entry per unique (sample, amplicon) combination that
passes all three filters.

**Header format:**

```
>SampleName    FwdTag1.RevTag1_FwdTag2.RevTag2    Count_rep1_Count_rep2
```

- Tags are joined by `.` within a replicate, and `_` across replicates.
- Counts are listed in replicate order.

**Example:**

```
>Sample1    tag1.tag2_tag3.tag4    24000_24000
AAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCC
>Sample2    tag5.tag6_tag7.tag8    5000_5000
AAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCC
```

---

## 8. Complex dataset walkthrough

This section demonstrates all sort and filter outcomes on a 100,000-read
synthetic dataset.

### 8.1 Dataset design

The dataset is generated by `test/generate_complex_dataset.py`.

**Primers (with IUPAC code):**

```
GCRTGC  AGTCAG
```

`R` is the IUPAC code for A or G, so the forward primer matches both
`GCATGC` and `GCGTGC`.  The dataset includes reads built with both variants
to demonstrate IUPAC-aware primer matching.

**Two samples, two PCR replicates each:**

| Sample | Fwd tag | Rev tag | Pool |
|--------|---------|---------|------|
| Sample1 | tag1 | tag2 | Pool1 (rep 1) |
| Sample1 | tag3 | tag4 | Pool2 (rep 2) |
| Sample2 | tag5 | tag6 | Pool1 (rep 1) |
| Sample2 | tag7 | tag8 | Pool2 (rep 2) |

**Amplicons in the dataset:**

| Name | Sequence | Length | Filter outcome |
|------|----------|--------|----------------|
| `AMP_PASSES_ALL` | `AAATTTGGGCCC` × 5 | 60 nt | PASSES (both reps, count ≥ 2, len ≥ 50) |
| `AMP_FAILS_PROP` | `GGGCCCAAATTT` × 5 | 60 nt | FAILS propPCRs=1.0 (only in Pool1) |
| `AMP_FAILS_COUNT` | `TTTGGGCCCAAA` × 5 | 60 nt | FAILS minOccurrence=2 (count=1 per rep) |
| `AMP_FAILS_LENGTH` | `GCTT` | 4 nt | FAILS minLength=50 (passes prop+count) |

**Sort outcomes represented:**

| Type | Description | Example |
|------|-------------|---------|
| C (forward) | Correct tag pair, forward orientation | tag1+tag2 + GCATGC |
| C (forward, IUPAC) | Correct tag pair, GCGTGC primer (R→G) | tag1+tag2 + GCGTGC |
| C (reverse) | Correct tag pair, reverse-oriented read | RC of tag1+tag2 read |
| B | Both tags in pool, wrong pair | tag1+tag6 in Pool1 |
| F | Forward tag in pool, reverse not | tag1+tag3 in Pool1 |
| R | Reverse tag in pool, forward not | tag3+tag2 in Pool1 |
| N | Neither tag in pool | tag3+tag4 in Pool1 |
| (none) | No primer — absent from tagInfo | random sequence |

### 8.2 Generate the dataset

```bash
python3 test/generate_complex_dataset.py
```

Output: `test/tutorial/` directory with all input files.

### 8.3 Run sort

**Python:**
```bash
mkdir -p test/tutorial/out
python3 begum/Begum.py sort \
  -p test/tutorial/primers.txt \
  -t test/tutorial/tags.txt \
  -s test/tutorial/samples.txt \
  -l test/tutorial/pools.txt \
  -d test/tutorial/out \
  -o tutorial
```

**Rust:**
```bash
mkdir -p test/tutorial/out
begum sort \
  -p test/tutorial/primers.txt \
  -t test/tutorial/tags.txt \
  -s test/tutorial/samples.txt \
  -l test/tutorial/pools.txt \
  -d test/tutorial/out \
  -o tutorial
```

### 8.4 Examine sort outputs

`test/tutorial/out/tutorial_Pool1.tagInfo` (key rows):

```
FTag    RTag    Seq                     Count   Type
tag1    tag2    AAATTTGGGCCC×5 (60nt)  24000   C   ← passes + 2000 reverse-orient reads
tag1    tag2    GCTT                    5000    C   ← will fail minLength
tag1    tag2    GGGCCCAAATTT×5 (60nt)  5000    C   ← will fail propPCRs (absent Pool2)
tag1    tag2    TTTGGGCCCAAA×5 (60nt)  1       C   ← will fail minOccurrence
tag5    tag6    AAATTTGGGCCC×5 (60nt)  5000    C   ← passes
tag5    tag6    GCTT                    1000    C   ← will fail minLength
tag5    tag6    TTTGGGCCCAAA×5 (60nt)  1       C   ← will fail minOccurrence
tag1    tag6    AAATTTGGGCCC×5 (60nt)  2000    B   ← wrong pair; ignored by filter
tag1    tag3    AAATTTGGGCCC×5 (60nt)  2000    F   ← fwd only; ignored by filter
tag3    tag2    AAATTTGGGCCC×5 (60nt)  2000    R   ← rev only; ignored by filter
tag3    tag4    AAATTTGGGCCC×5 (60nt)  1000    N   ← neither; ignored by filter
```

Pool2 is similar but has tag3/tag4 (Sample1 rep2) and tag7/tag8 (Sample2 rep2).
`AMP_FAILS_PROP` is absent from Pool2, confirming it belongs only to Pool1.

`test/tutorial/out/tutorial_Pool1.summaryCounts`:

```
Tag1    Tag2    UniqSeqs    TotalSeqs   Type
Correct combination of tags used in pool
----------------------------------------
tag1    tag2    4           34001       C
tag5    tag6    3           6001        C
Both tags used in pool, but not this combination
------------------------------------------------
tag1    tag6    1           2000        B
Only forward tag used in pool, reverse not found
------------------------------------------------
tag1    tag3    1           2000        F
Only reverse tag used in pool, forward not found
------------------------------------------------
tag3    tag2    1           2000        R
Neither tag used in this pool
-----------------------------
tag3    tag4    1           1000        N
```

### 8.5 Run filter

**Python:**
```bash
python3 begum/Begum.py filter \
  -i test/tutorial/out/tutorial \
  -s test/tutorial/samples.txt \
  -p 1.0 \
  -m 2 \
  -l 50 \
  -d test/tutorial/out \
  -o filtered
```

**Rust:**
```bash
begum filter \
  -i test/tutorial/out/tutorial \
  -s test/tutorial/samples.txt \
  -p 1.0 \
  -m 2 \
  -l 50 \
  -d test/tutorial/out \
  -o filtered
```

### 8.6 Examine filter output

`test/tutorial/out/filtered.fna`:

```
>Sample1    tag1.tag2_tag3.tag4    24000_24000
AAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCC
>Sample2    tag5.tag6_tag7.tag8    5000_5000
AAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCC
```

Only `AMP_PASSES_ALL` survives for both samples (header fields are tab-separated).  The Rust version may list replicates in a different order in the header (e.g. `tag7.tag8_tag5.tag6` instead of `tag5.tag6_tag7.tag8`) due to HashMap iteration order, but sequence content and counts are identical.

The other amplicons are excluded for the following reasons:

| Amplicon | Sample | Count Pool1 | Count Pool2 | Reps passed | Outcome |
|----------|--------|-------------|-------------|-------------|---------|
| `AMP_PASSES_ALL` | Sample1 | 24000 | 24000 | 2/2 = 1.0 | ✓ PASSES |
| `AMP_PASSES_ALL` | Sample2 | 5000 | 5000 | 2/2 = 1.0 | ✓ PASSES |
| `AMP_FAILS_PROP` | Sample1 | 5000 | 0 | 1/2 = 0.5 < 1.0 | ✗ propPCRs |
| `AMP_FAILS_COUNT` | Sample1 | 1 (<2) | 1 (<2) | 0/2 = 0.0 | ✗ minOccurrence |
| `AMP_FAILS_COUNT` | Sample2 | 1 (<2) | 1 (<2) | 0/2 = 0.0 | ✗ minOccurrence |
| `AMP_FAILS_LENGTH` | Sample1 | 5000 | 5000 | 2/2 = 1.0 | ✗ minLength (4 < 50) |
| `AMP_FAILS_LENGTH` | Sample2 | 1000 | 1000 | 2/2 = 1.0 | ✗ minLength (4 < 50) |

The count of 24000 for Sample1 Pool1 reflects 22,000 standard reads +
2,000 IUPAC-matched reads (GCGTGC primer) + 2,000 reverse-orientation reads
— all mapping to the same amplicon and tag pair.
