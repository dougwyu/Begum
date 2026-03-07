use begum::{filter, sort, FilterArgs, SortArgs};
use std::{fs, io::Write, path::PathBuf};
use tempfile::tempdir;

fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            _ => b'N',
        })
        .collect()
}

fn make_read(fwd_tag: &[u8], rev_tag: &[u8], fwd_primer: &[u8], rev_primer: &[u8], amplicon: &[u8]) -> Vec<u8> {
    let rc_rev_primer = rc(rev_primer);
    let rc_rev_tag = rc(rev_tag);
    [fwd_tag, fwd_primer, amplicon, &rc_rev_primer, &rc_rev_tag].concat()
}

fn write_fastq(path: &PathBuf, reads: &[Vec<u8>]) {
    let mut f = fs::File::create(path).unwrap();
    for (i, seq) in reads.iter().enumerate() {
        let qual = "I".repeat(seq.len());
        writeln!(
            f,
            "@read{}\n{}\n+\n{}",
            i + 1,
            std::str::from_utf8(seq).unwrap(),
            qual
        )
        .unwrap();
    }
}

#[test]
fn test_sort_and_filter_end_to_end() {
    let dir = tempdir().unwrap();
    let d = dir.path();

    let fwd_primer = b"GCATGC";
    let rev_primer = b"AGTCAG";
    let tag1 = b"AACCGGT";
    let tag2 = b"TTGGCCA";
    let tag3 = b"CCGGAAT";
    let tag4 = b"GGCCTTA";
    let amp_shared = b"AAATTTGGGCCC";
    let amp_pool1_only = b"GGGCCCAAATTT";

    // Pool1: 2x shared, 1x pool1-only
    let pool1_reads = vec![
        make_read(tag1, tag2, fwd_primer, rev_primer, amp_shared),
        make_read(tag1, tag2, fwd_primer, rev_primer, amp_shared),
        make_read(tag1, tag2, fwd_primer, rev_primer, amp_pool1_only),
    ];
    write_fastq(&d.join("reads_pool1.fastq"), &pool1_reads);

    // Pool2: 1x shared only
    let pool2_reads = vec![make_read(tag3, tag4, fwd_primer, rev_primer, amp_shared)];
    write_fastq(&d.join("reads_pool2.fastq"), &pool2_reads);

    // Write input files
    fs::write(d.join("primers.txt"), "GCATGC AGTCAG\n").unwrap();
    fs::write(
        d.join("tags.txt"),
        "tag1 AACCGGT\ntag2 TTGGCCA\ntag3 CCGGAAT\ntag4 GGCCTTA\n",
    )
    .unwrap();
    fs::write(
        d.join("samples.txt"),
        "Sample1 tag1 tag2 Pool1\nSample1 tag3 tag4 Pool2\n",
    )
    .unwrap();
    fs::write(
        d.join("pools.txt"),
        format!(
            "Pool1 {}\nPool2 {}\n",
            d.join("reads_pool1.fastq").display(),
            d.join("reads_pool2.fastq").display()
        ),
    )
    .unwrap();

    // ── Run sort ──────────────────────────────────────────────────────────
    sort::run(SortArgs {
        primers: Some(d.join("primers.txt").to_str().unwrap().to_string()),
        fwd_primer: None,
        rev_primer: None,
        tags: d.join("tags.txt").to_str().unwrap().to_string(),
        sample_info: d.join("samples.txt").to_str().unwrap().to_string(),
        pool: d.join("pools.txt").to_str().unwrap().to_string(),
        allow_multiple_primers: false,
        primer_mismatches: 0,
        tag_mismatches: 0,
        output_directory: d.to_str().unwrap().to_string(),
        output_prefix: "test".to_string(),
    })
    .unwrap();

    // Verify Pool1 tagInfo
    let p1 = fs::read_to_string(d.join("test_Pool1.tagInfo")).unwrap();
    let p1_rows: Vec<Vec<&str>> = p1.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    let shared_row = p1_rows.iter().find(|r| r.get(2) == Some(&"AAATTTGGGCCC")).unwrap();
    assert_eq!(shared_row[3], "2", "shared amplicon count=2 in Pool1");
    assert_eq!(shared_row[4], "C", "shared amplicon type=C in Pool1");
    let rare_row = p1_rows.iter().find(|r| r.get(2) == Some(&"GGGCCCAAATTT")).unwrap();
    assert_eq!(rare_row[3], "1", "pool1-only count=1");

    // Verify Pool2 tagInfo
    let p2 = fs::read_to_string(d.join("test_Pool2.tagInfo")).unwrap();
    let p2_rows: Vec<Vec<&str>> = p2.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    assert_eq!(p2_rows.len(), 1, "Pool2 has exactly 1 amplicon");
    assert_eq!(p2_rows[0][2], "AAATTTGGGCCC", "shared amplicon in Pool2");
    assert_eq!(p2_rows[0][3], "1", "count=1 in Pool2");

    // ── Run filter ────────────────────────────────────────────────────────
    filter::run_filter(
        d.join("test").to_str().unwrap(),
        d.join("samples.txt").to_str().unwrap(),
        1.0, 1, 5,
        d.to_str().unwrap(),
        "filtered",
    )
    .unwrap();

    let fna = fs::read_to_string(d.join("filtered.fna")).unwrap();
    let seqs: Vec<&str> = fna.lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(seqs.len(), 1, "only 1 sequence survives propPCRs=1.0");
    assert!(
        seqs.contains(&"AAATTTGGGCCC"),
        "shared amplicon must be in output"
    );
    assert!(
        !fna.contains("GGGCCCAAATTT"),
        "pool1-only amplicon must be filtered out"
    );
}
