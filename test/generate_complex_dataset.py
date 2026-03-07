#!/usr/bin/env python3
"""generate_complex_dataset.py

Generate 100,000 synthetic FASTQ reads across 2 pools that demonstrate every
sort and filter outcome.  Outputs all input files needed to run Begum in the
test/tutorial/ directory.

Usage:
    python3 test/generate_complex_dataset.py

Outputs (relative to repo root):
    test/tutorial/reads_pool1.fastq   (50,000 reads)
    test/tutorial/reads_pool2.fastq   (50,000 reads)
    test/tutorial/primers.txt
    test/tutorial/tags.txt
    test/tutorial/samples.txt
    test/tutorial/pools.txt
"""

import gzip
import random
from pathlib import Path

SEED = 42
random.seed(SEED)

# ---------------------------------------------------------------------------
# Experiment design
# ---------------------------------------------------------------------------
#
#  Primers (IUPAC R = A or G):
#    Forward:  GCRTGC   matches both GCATGC and GCGTGC
#    Reverse:  AGTCAG
#
#  Samples and PCR replicates:
#    Sample1  tag1 (fwd) + tag2 (rev)  Pool1   → replicate 1
#    Sample1  tag3 (fwd) + tag4 (rev)  Pool2   → replicate 2
#    Sample2  tag5 (fwd) + tag6 (rev)  Pool1   → replicate 1
#    Sample2  tag7 (fwd) + tag8 (rev)  Pool2   → replicate 2
#
#  Tags are 7 nt, pairwise Hamming distance ≥ 3.
#
#  Amplicons (filter run with -p 1.0 -m 2 -l 50):
#    AMP_PASSES_ALL   60 nt  in both reps, count>>2, len>>50  → PASSES all
#    AMP_FAILS_PROP   60 nt  only in Pool1                    → FAILS propPCRs
#    AMP_FAILS_COUNT  60 nt  count=1 per rep                  → FAILS minOccurrence
#    AMP_FAILS_LENGTH  4 nt  count>>2 in both reps            → FAILS minLength
#
#  Sort types in tagInfo (Type column):
#    C  correct tag combination for that pool
#    B  both tags present in pool but wrong pair
#    F  only the forward tag found in pool
#    R  only the reverse tag found in pool
#    N  neither tag found in pool
#    (reads with no primer match do not appear in tagInfo at all)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FWD_PRIMER_IUPAC = "GCRTGC"
REV_PRIMER_STR   = "AGTCAG"

FWD_PRIMER_A = b"GCATGC"   # IUPAC R resolved to A
FWD_PRIMER_G = b"GCGTGC"   # IUPAC R resolved to G — demonstrates IUPAC matching
RC_REV_PRIMER = b"CTGACT"  # reverse complement of AGTCAG

TAGS = {
    "tag1": b"AACCGGT",   # Sample1 | Pool1 | fwd
    "tag2": b"TTGGCCA",   # Sample1 | Pool1 | rev
    "tag3": b"CCGGAAT",   # Sample1 | Pool2 | fwd
    "tag4": b"GGCCTTA",   # Sample1 | Pool2 | rev
    "tag5": b"AATCCGG",   # Sample2 | Pool1 | fwd
    "tag6": b"TTAAGGC",   # Sample2 | Pool1 | rev
    "tag7": b"GGTAACC",   # Sample2 | Pool2 | fwd
    "tag8": b"CCAATTG",   # Sample2 | Pool2 | rev
}

# Amplicons
AMP_PASSES_ALL   = b"AAATTTGGGCCC" * 5   # 60 nt
AMP_FAILS_PROP   = b"GGGCCCAAATTT" * 5   # 60 nt, only in Pool1
AMP_FAILS_COUNT  = b"TTTGGGCCCAAA" * 5   # 60 nt, count=1 per rep
AMP_FAILS_LENGTH = b"GCTT"               # 4 nt

# Used as amplicon for B/F/R/N reads (amplicon content irrelevant to sort type)
_DUMMY_AMP = AMP_PASSES_ALL

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

_COMP = {65: 84, 84: 65, 71: 67, 67: 71}  # A↔T, G↔C (ASCII values)

def rc(seq: bytes) -> bytes:
    """Reverse complement of a DNA byte string."""
    return bytes(_COMP.get(b, 78) for b in reversed(seq))


def make_fwd_read(fwd_tag: bytes, fwd_primer: bytes,
                  amplicon: bytes, rev_tag: bytes) -> bytes:
    """Build a forward-orientation read.

    Structure: [fwd_tag][fwd_primer][amplicon][rc(rev_primer)][rc(rev_tag)]
    """
    return fwd_tag + fwd_primer + amplicon + RC_REV_PRIMER + rc(rev_tag)


def make_rev_read(fwd_tag: bytes, fwd_primer: bytes,
                  amplicon: bytes, rev_tag: bytes) -> bytes:
    """Build a reverse-orientation (type-4) read.

    This is the reverse complement of the forward read.  The sort step
    reverse-complements the read and finds the same tag/primer structure,
    ultimately extracting the same tag combo and amplicon.
    """
    return rc(make_fwd_read(fwd_tag, fwd_primer, amplicon, rev_tag))


def random_noprim_seq(length: int = 150) -> bytes:
    """Random sequence guaranteed not to contain any primer or its RC."""
    avoid = {b"GCATGC", b"GCGTGC", b"CTGACT", b"AGTCAG"}
    while True:
        seq = bytes(random.choice(b"ACGT") for _ in range(length))
        if not any(p in seq for p in avoid):
            return seq


def write_fastq(path: Path, reads: list) -> None:
    """Write a list of byte sequences as a gzip-compressed FASTQ file."""
    with gzip.open(path, "wt") as f:
        for i, seq in enumerate(reads, 1):
            qual = "I" * len(seq)
            f.write(f"@read{i}\n{seq.decode()}\n+\n{qual}\n")


# ---------------------------------------------------------------------------
# Build Pool1 reads
# ---------------------------------------------------------------------------
# Pool1 tags: tag1, tag2 (Sample1 rep1) and tag5, tag6 (Sample2 rep1)
# Valid pairs in Pool1: (tag1, tag2) and (tag5, tag6)
#
# B-type tag combo for Pool1: (tag1, tag6) — both in pool but wrong pair
# F-type tag combo for Pool1: (tag1, tag3) — tag1 in pool, tag3 NOT in pool
# R-type tag combo for Pool1: (tag3, tag2) — tag3 NOT in pool, tag2 in pool
# N-type tag combo for Pool1: (tag3, tag4) — neither in pool

pool1_reads = []

# --- Type C, correct combinations -------------------------------------------

# Sample1 rep1, amp_passes_all, primer variant A (standard)
pool1_reads += [make_fwd_read(TAGS["tag1"], FWD_PRIMER_A, AMP_PASSES_ALL, TAGS["tag2"])
                for _ in range(22_000)]

# Sample1 rep1, amp_passes_all, primer variant G  →  demonstrates IUPAC primer
pool1_reads += [make_fwd_read(TAGS["tag1"], FWD_PRIMER_G, AMP_PASSES_ALL, TAGS["tag2"])
                for _ in range(2_000)]

# Sample1 rep1, amp_fails_prop (amplicon absent from Pool2 → fails propPCRs=1.0)
pool1_reads += [make_fwd_read(TAGS["tag1"], FWD_PRIMER_A, AMP_FAILS_PROP, TAGS["tag2"])
                for _ in range(5_000)]

# Sample1 rep1, amp_fails_count (count=1 per rep → fails minOccurrence=2)
pool1_reads.append(make_fwd_read(TAGS["tag1"], FWD_PRIMER_A, AMP_FAILS_COUNT, TAGS["tag2"]))

# Sample1 rep1, amp_fails_length (4 nt amplicon → fails minLength=50)
pool1_reads += [make_fwd_read(TAGS["tag1"], FWD_PRIMER_A, AMP_FAILS_LENGTH, TAGS["tag2"])
                for _ in range(5_000)]

# Sample2 rep1, amp_passes_all
pool1_reads += [make_fwd_read(TAGS["tag5"], FWD_PRIMER_A, AMP_PASSES_ALL, TAGS["tag6"])
                for _ in range(5_000)]

# Sample2 rep1, amp_fails_count
pool1_reads.append(make_fwd_read(TAGS["tag5"], FWD_PRIMER_A, AMP_FAILS_COUNT, TAGS["tag6"]))

# Sample2 rep1, amp_fails_length
pool1_reads += [make_fwd_read(TAGS["tag5"], FWD_PRIMER_A, AMP_FAILS_LENGTH, TAGS["tag6"])
                for _ in range(1_000)]

# --- Type C, reverse orientation (type-4) ------------------------------------
# Same sample/amplicon but read is RC of normal read.  Sort extracts the same
# tag combo and amplicon from the reverse complement of the input read.

pool1_reads += [make_rev_read(TAGS["tag1"], FWD_PRIMER_A, AMP_PASSES_ALL, TAGS["tag2"])
                for _ in range(2_000)]

# --- Type B: both tags in pool, but not the expected pair --------------------
pool1_reads += [make_fwd_read(TAGS["tag1"], FWD_PRIMER_A, _DUMMY_AMP, TAGS["tag6"])
                for _ in range(2_000)]

# --- Type F: fwd tag in pool, rev tag not in pool ----------------------------
pool1_reads += [make_fwd_read(TAGS["tag1"], FWD_PRIMER_A, _DUMMY_AMP, TAGS["tag3"])
                for _ in range(2_000)]

# --- Type R: fwd tag not in pool, rev tag in pool ----------------------------
pool1_reads += [make_fwd_read(TAGS["tag3"], FWD_PRIMER_A, _DUMMY_AMP, TAGS["tag2"])
                for _ in range(2_000)]

# --- Type N: neither tag in pool ---------------------------------------------
pool1_reads += [make_fwd_read(TAGS["tag3"], FWD_PRIMER_A, _DUMMY_AMP, TAGS["tag4"])
                for _ in range(1_000)]

# --- No-primer reads (don't appear in tagInfo at all) -----------------------
pool1_reads += [random_noprim_seq() for _ in range(998)]

assert len(pool1_reads) == 50_000, f"Pool1 has {len(pool1_reads)} reads, expected 50,000"

random.shuffle(pool1_reads)

# ---------------------------------------------------------------------------
# Build Pool2 reads
# ---------------------------------------------------------------------------
# Pool2 tags: tag3, tag4 (Sample1 rep2) and tag7, tag8 (Sample2 rep2)
# Valid pairs in Pool2: (tag3, tag4) and (tag7, tag8)
#
# B-type tag combo for Pool2: (tag3, tag8)
# F-type tag combo for Pool2: (tag3, tag1) — tag1 NOT in pool2
# R-type tag combo for Pool2: (tag1, tag4) — tag1 NOT in pool2
# N-type tag combo for Pool2: (tag1, tag2) — neither in pool2

pool2_reads = []

# --- Type C, correct combinations -------------------------------------------

# Sample1 rep2, amp_passes_all, primer variant A
pool2_reads += [make_fwd_read(TAGS["tag3"], FWD_PRIMER_A, AMP_PASSES_ALL, TAGS["tag4"])
                for _ in range(22_000)]

# Sample1 rep2, amp_passes_all, primer variant G (IUPAC demo)
pool2_reads += [make_fwd_read(TAGS["tag3"], FWD_PRIMER_G, AMP_PASSES_ALL, TAGS["tag4"])
                for _ in range(2_000)]

# Sample1 rep2, amp_fails_count
pool2_reads.append(make_fwd_read(TAGS["tag3"], FWD_PRIMER_A, AMP_FAILS_COUNT, TAGS["tag4"]))

# Sample1 rep2, amp_fails_length
pool2_reads += [make_fwd_read(TAGS["tag3"], FWD_PRIMER_A, AMP_FAILS_LENGTH, TAGS["tag4"])
                for _ in range(5_000)]

# Sample2 rep2, amp_passes_all
pool2_reads += [make_fwd_read(TAGS["tag7"], FWD_PRIMER_A, AMP_PASSES_ALL, TAGS["tag8"])
                for _ in range(5_000)]

# Sample2 rep2, amp_fails_count
pool2_reads.append(make_fwd_read(TAGS["tag7"], FWD_PRIMER_A, AMP_FAILS_COUNT, TAGS["tag8"]))

# Sample2 rep2, amp_fails_length
pool2_reads += [make_fwd_read(TAGS["tag7"], FWD_PRIMER_A, AMP_FAILS_LENGTH, TAGS["tag8"])
                for _ in range(1_000)]

# --- Type C, reverse orientation (type-4) ------------------------------------
pool2_reads += [make_rev_read(TAGS["tag3"], FWD_PRIMER_A, AMP_PASSES_ALL, TAGS["tag4"])
                for _ in range(2_000)]

# --- Type B -------------------------------------------------------------------
pool2_reads += [make_fwd_read(TAGS["tag3"], FWD_PRIMER_A, _DUMMY_AMP, TAGS["tag8"])
                for _ in range(3_000)]

# --- Type F -------------------------------------------------------------------
pool2_reads += [make_fwd_read(TAGS["tag3"], FWD_PRIMER_A, _DUMMY_AMP, TAGS["tag1"])
                for _ in range(3_000)]

# --- Type R -------------------------------------------------------------------
pool2_reads += [make_fwd_read(TAGS["tag1"], FWD_PRIMER_A, _DUMMY_AMP, TAGS["tag4"])
                for _ in range(3_000)]

# --- Type N -------------------------------------------------------------------
pool2_reads += [make_fwd_read(TAGS["tag1"], FWD_PRIMER_A, _DUMMY_AMP, TAGS["tag2"])
                for _ in range(2_000)]

# --- No-primer reads ----------------------------------------------------------
pool2_reads += [random_noprim_seq() for _ in range(1_998)]

assert len(pool2_reads) == 50_000, f"Pool2 has {len(pool2_reads)} reads, expected 50,000"

random.shuffle(pool2_reads)

# ---------------------------------------------------------------------------
# Write output files
# ---------------------------------------------------------------------------

OUT = Path(__file__).parent / "tutorial"
OUT.mkdir(exist_ok=True)

write_fastq(OUT / "reads_pool1.fastq.gz", pool1_reads)
write_fastq(OUT / "reads_pool2.fastq.gz", pool2_reads)
print(f"Wrote {len(pool1_reads):,} reads to {OUT / 'reads_pool1.fastq.gz'}")
print(f"Wrote {len(pool2_reads):,} reads to {OUT / 'reads_pool2.fastq.gz'}")

# primers.txt — single line: fwd_primer rev_primer
(OUT / "primers.txt").write_text(f"{FWD_PRIMER_IUPAC} {REV_PRIMER_STR}\n")

# tags.txt — one tag per line: name  sequence
tags_content = "\n".join(f"{name} {seq.decode()}" for name, seq in TAGS.items()) + "\n"
(OUT / "tags.txt").write_text(tags_content)

# samples.txt — one replicate per line: SampleName FwdTag RevTag PoolName
(OUT / "samples.txt").write_text(
    "Sample1 tag1 tag2 Pool1\n"
    "Sample1 tag3 tag4 Pool2\n"
    "Sample2 tag5 tag6 Pool1\n"
    "Sample2 tag7 tag8 Pool2\n"
)

# pools.txt — one pool per line: PoolName FASTQ_path
pool1_path = (OUT / "reads_pool1.fastq.gz").resolve()
pool2_path = (OUT / "reads_pool2.fastq.gz").resolve()
(OUT / "pools.txt").write_text(
    f"Pool1 {pool1_path}\n"
    f"Pool2 {pool2_path}\n"
)

print(f"Wrote input files to {OUT}/")
print()
print("Expected filter results (with -p 1.0 -m 2 -l 50):")
print("  Sample1: AMP_PASSES_ALL (60 nt) — PASSES all filters")
print("  Sample2: AMP_PASSES_ALL (60 nt) — PASSES all filters")
print("  Sample1: AMP_FAILS_PROP (60 nt) — FAILS propPCRs (only in Pool1)")
print("  Sample1: AMP_FAILS_COUNT (60 nt) — FAILS minOccurrence (count=1 per rep)")
print("  Sample2: AMP_FAILS_COUNT (60 nt) — FAILS minOccurrence (count=1 per rep)")
print("  Sample1: AMP_FAILS_LENGTH (4 nt) — FAILS minLength (4 < 50)")
print("  Sample2: AMP_FAILS_LENGTH (4 nt) — FAILS minLength (4 < 50)")
