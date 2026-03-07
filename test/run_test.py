#!/usr/bin/env python3
"""End-to-end test for Begum sort and filter."""

import gzip
import sys
import subprocess
from pathlib import Path

TEST_DIR = Path(__file__).parent.resolve()
SRC_DIR = TEST_DIR.parent / "begum"
OUT_DIR = TEST_DIR / "out"

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  [{PASS}] {name}")
    else:
        print(f"  [{FAIL}] {name}" + (f": {detail}" if detail else ""))
        failures.append(name)


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def rc(seq):
    return seq.translate(_COMP)[::-1]

def make_read(fwd_tag, rev_tag, fwd_primer, rev_primer, amplicon):
    return fwd_tag + fwd_primer + amplicon + rc(rev_primer) + rc(rev_tag)


def write_fastq(path, reads):
    with gzip.open(path, "wt") as f:
        for i, seq in enumerate(reads, 1):
            f.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def run_begum(args):
    cmd = [sys.executable, str(SRC_DIR / "Begum.py")] + args
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SRC_DIR))
    return result


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
OUT_DIR.mkdir(exist_ok=True)

FWD_PRIMER = "GCATGC"
REV_PRIMER = "AGTCAG"

# Pool1 tags → Sample1 replicate 1
TAG1, TAG2 = "AACCGGT", "TTGGCCA"
# Pool2 tags → Sample1 replicate 2
TAG3, TAG4 = "CCGGAAT", "GGCCTTA"

AMP_SHARED    = "AAATTTGGGCCC"   # in both pools → passes propPCRs=1.0
AMP_POOL1_ONLY = "GGGCCCAAATTT"  # only in Pool1 → filtered by propPCRs=1.0

# Pool1: 2x shared amplicon, 1x pool1-only amplicon
pool1_reads = [
    make_read(TAG1, TAG2, FWD_PRIMER, REV_PRIMER, AMP_SHARED),
    make_read(TAG1, TAG2, FWD_PRIMER, REV_PRIMER, AMP_SHARED),
    make_read(TAG1, TAG2, FWD_PRIMER, REV_PRIMER, AMP_POOL1_ONLY),
]
write_fastq(TEST_DIR / "reads_pool1.fastq.gz", pool1_reads)

# Pool2: 1x shared amplicon only
pool2_reads = [
    make_read(TAG3, TAG4, FWD_PRIMER, REV_PRIMER, AMP_SHARED),
]
write_fastq(TEST_DIR / "reads_pool2.fastq.gz", pool2_reads)

(TEST_DIR / "primers.txt").write_text(f"{FWD_PRIMER} {REV_PRIMER}\n")
(TEST_DIR / "tags.txt").write_text(
    f"tag1 {TAG1}\ntag2 {TAG2}\ntag3 {TAG3}\ntag4 {TAG4}\n"
)
(TEST_DIR / "samples.txt").write_text(
    "Sample1 tag1 tag2 Pool1\n"
    "Sample1 tag3 tag4 Pool2\n"
)
(TEST_DIR / "pools.txt").write_text(
    f"Pool1 {TEST_DIR}/reads_pool1.fastq.gz\n"
    f"Pool2 {TEST_DIR}/reads_pool2.fastq.gz\n"
)

# ---------------------------------------------------------------------------
# Test: sort
# ---------------------------------------------------------------------------
print("\n=== sort ===")
result = run_begum([
    "sort",
    "-p", str(TEST_DIR / "primers.txt"),
    "-t", str(TEST_DIR / "tags.txt"),
    "-s", str(TEST_DIR / "samples.txt"),
    "-l", str(TEST_DIR / "pools.txt"),
    "-d", str(OUT_DIR),
    "-o", "test",
])
check("sort exits 0", result.returncode == 0, result.stderr)

# Parse Pool1 tagInfo
tag_info_pool1 = OUT_DIR / "test_Pool1.tagInfo"
check("Pool1 tagInfo exists", tag_info_pool1.exists())

if tag_info_pool1.exists():
    lines = tag_info_pool1.read_text().splitlines()
    data = {row[2]: int(row[3]) for row in (l.split("\t") for l in lines[1:]) if row[4] == "C"}
    check("Pool1: shared amplicon count=2",  data.get(AMP_SHARED) == 2,
          f"got {data.get(AMP_SHARED)}")
    check("Pool1: pool1-only amplicon count=1", data.get(AMP_POOL1_ONLY) == 1,
          f"got {data.get(AMP_POOL1_ONLY)}")
    check("Pool1: only 2 amplicons in tagInfo", len(data) == 2, f"got {len(data)}")

# Parse Pool2 tagInfo
tag_info_pool2 = OUT_DIR / "test_Pool2.tagInfo"
check("Pool2 tagInfo exists", tag_info_pool2.exists())

if tag_info_pool2.exists():
    lines = tag_info_pool2.read_text().splitlines()
    data2 = {row[2]: int(row[3]) for row in (l.split("\t") for l in lines[1:]) if row[4] == "C"}
    check("Pool2: shared amplicon count=1", data2.get(AMP_SHARED) == 1,
          f"got {data2.get(AMP_SHARED)}")
    check("Pool2: pool1-only amplicon absent", AMP_POOL1_ONLY not in data2)
    check("Pool2: only 1 amplicon in tagInfo", len(data2) == 1, f"got {len(data2)}")

# ---------------------------------------------------------------------------
# Test: filter  (propPCRs=1.0, minOccurence=1, minLength=5)
# ---------------------------------------------------------------------------
print("\n=== filter ===")
result = run_begum([
    "filter",
    "-i", str(OUT_DIR / "test"),
    "-s", str(TEST_DIR / "samples.txt"),
    "-p", "1.0",
    "-m", "1",
    "-l", "5",
    "-d", str(OUT_DIR),
    "-o", "filtered",
])
check("filter exits 0", result.returncode == 0, result.stderr)

fna = OUT_DIR / "filtered.fna"
check("filtered.fna exists", fna.exists())

if fna.exists():
    fna_text = fna.read_text()
    seqs = [l for l in fna_text.splitlines() if not l.startswith(">")]
    check("shared amplicon in output",    AMP_SHARED in seqs)
    check("pool1-only amplicon filtered", AMP_POOL1_ONLY not in seqs,
          "propPCRs=1.0 should have removed it")
    check("exactly 1 sequence in output", len(seqs) == 1, f"got {len(seqs)}")

# ---------------------------------------------------------------------------
# Setup: paired-end
# ---------------------------------------------------------------------------
PE_OUT_DIR = TEST_DIR / "out_pe"
PE_OUT_DIR.mkdir(exist_ok=True)

# Paired-end read structure:
#   Read1 = fwd_tag + fwd_primer + amplicon_R1
#   Read2 = rev_tag + rev_primer + amplicon_R2
# (no RC needed: the two reads are sequenced from opposite ends)
AMP_SHARED_R1   = "AAATTTGGG"
AMP_SHARED_R2   = "CCCTTTAGG"
AMP_POOL1_R1    = "GGGCCCAAA"
AMP_POOL1_R2    = "TTTGGGCCC"

pool1_r1 = [
    TAG1 + FWD_PRIMER + AMP_SHARED_R1,
    TAG1 + FWD_PRIMER + AMP_SHARED_R1,
    TAG1 + FWD_PRIMER + AMP_POOL1_R1,
]
pool1_r2 = [
    TAG2 + REV_PRIMER + AMP_SHARED_R2,
    TAG2 + REV_PRIMER + AMP_SHARED_R2,
    TAG2 + REV_PRIMER + AMP_POOL1_R2,
]
pool2_r1 = [TAG3 + FWD_PRIMER + AMP_SHARED_R1]
pool2_r2 = [TAG4 + REV_PRIMER + AMP_SHARED_R2]

write_fastq(TEST_DIR / "pe_pool1_R1.fastq.gz", pool1_r1)
write_fastq(TEST_DIR / "pe_pool1_R2.fastq.gz", pool1_r2)
write_fastq(TEST_DIR / "pe_pool2_R1.fastq.gz", pool2_r1)
write_fastq(TEST_DIR / "pe_pool2_R2.fastq.gz", pool2_r2)

(TEST_DIR / "pe_primers.txt").write_text(f"{FWD_PRIMER} {REV_PRIMER}\n")
(TEST_DIR / "pe_tags.txt").write_text(
    f"tag1 {TAG1}\ntag2 {TAG2}\ntag3 {TAG3}\ntag4 {TAG4}\n"
)
(TEST_DIR / "pe_samples.txt").write_text(
    "Sample1 tag1 tag2 Pool1\n"
    "Sample1 tag3 tag4 Pool2\n"
)
(TEST_DIR / "pe_pools.txt").write_text(
    f"Pool1 {TEST_DIR}/pe_pool1_R1.fastq.gz {TEST_DIR}/pe_pool1_R2.fastq.gz\n"
    f"Pool2 {TEST_DIR}/pe_pool2_R1.fastq.gz {TEST_DIR}/pe_pool2_R2.fastq.gz\n"
)

# ---------------------------------------------------------------------------
# Test: sort (paired-end)
# ---------------------------------------------------------------------------
print("\n=== sort (paired-end) ===")
result = run_begum([
    "sort",
    "-p", str(TEST_DIR / "pe_primers.txt"),
    "-t", str(TEST_DIR / "pe_tags.txt"),
    "-s", str(TEST_DIR / "pe_samples.txt"),
    "-l", str(TEST_DIR / "pe_pools.txt"),
    "-d", str(PE_OUT_DIR),
    "-o", "test_pe",
])
check("PE sort exits 0", result.returncode == 0, result.stderr)

pe_tag_info = PE_OUT_DIR / "test_pe_Pool1.tagInfo"
check("PE Pool1 tagInfo exists", pe_tag_info.exists())

if pe_tag_info.exists():
    lines = pe_tag_info.read_text().splitlines()
    header_cols = lines[0].split("\t")
    check("PE tagInfo has 6 columns", len(header_cols) == 6,
          f"got {len(header_cols)}")
    # 6-column format: FTag RTag FSeq RSeq Count Type
    data_pe = {
        (row[2], row[3]): int(row[4])
        for row in (l.split("\t") for l in lines[1:])
        if row[5] == "C"
    }
    check("PE Pool1: shared amplicon count=2",
          data_pe.get((AMP_SHARED_R1, AMP_SHARED_R2)) == 2,
          f"got {data_pe.get((AMP_SHARED_R1, AMP_SHARED_R2))}")
    check("PE Pool1: pool1-only amplicon count=1",
          data_pe.get((AMP_POOL1_R1, AMP_POOL1_R2)) == 1,
          f"got {data_pe.get((AMP_POOL1_R1, AMP_POOL1_R2))}")

# ---------------------------------------------------------------------------
# Test: filter (paired-end)  (propPCRs=1.0, minOccurrence=1, minLength=5)
# ---------------------------------------------------------------------------
print("\n=== filter (paired-end) ===")
result = run_begum([
    "filter",
    "-i", str(PE_OUT_DIR / "test_pe"),
    "-s", str(TEST_DIR / "pe_samples.txt"),
    "-p", "1.0",
    "-m", "1",
    "-l", "5",
    "-d", str(PE_OUT_DIR),
    "-o", "filtered_pe",
])
check("PE filter exits 0", result.returncode == 0, result.stderr)

pe_fna = PE_OUT_DIR / "filtered_pe.fna"
check("PE filtered.fna exists", pe_fna.exists())

if pe_fna.exists():
    fna_text = pe_fna.read_text()
    seqs = [l for l in fna_text.splitlines() if not l.startswith(">")]
    headers = [l for l in fna_text.splitlines() if l.startswith(">")]
    check("PE shared R1 in output",      AMP_SHARED_R1 in seqs)
    check("PE shared R2 in output",      AMP_SHARED_R2 in seqs)
    check("PE pool1-only R1 filtered",   AMP_POOL1_R1 not in seqs)
    check("PE pool1-only R2 filtered",   AMP_POOL1_R2 not in seqs)
    check("PE exactly 2 seqs in output", len(seqs) == 2, f"got {len(seqs)}")
    check("PE headers have read1/read2",
          any("read1" in h for h in headers) and any("read2" in h for h in headers))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print()
if failures:
    print(f"FAILED ({len(failures)} checks): {', '.join(failures)}")
    sys.exit(1)
else:
    print(f"All checks passed.")
