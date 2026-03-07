use crate::args::FilterArgs;
use anyhow::{Context, Result};
use std::{collections::HashMap, fs, io::Write};

// ── Types ──────────────────────────────────────────────────────────────────

/// sample → replicate → "fwd_tag.rev_tag" label
type RepInfo = HashMap<String, HashMap<usize, String>>;

/// sample → seq → replicate → (tag_pair_label, count)
type HapsInfo = HashMap<String, HashMap<String, HashMap<usize, (String, u32)>>>;

// ── Public entry points ────────────────────────────────────────────────────

pub fn run(args: FilterArgs) -> Result<()> {
    run_filter(
        &args.input_prefix,
        &args.sample_info,
        args.prop_pcrs,
        args.min_occurrence,
        args.min_length,
        &args.output_directory,
        &args.output_prefix,
    )
}

pub fn run_filter(
    input_prefix: &str,
    sample_info_path: &str,
    prop_pcrs: f64,
    min_occurrence: usize,
    min_length: usize,
    output_directory: &str,
    output_prefix: &str,
) -> Result<()> {
    // Re-use sort's sample info reader
    let sample_info = crate::sort::read_sample_info_file(sample_info_path)?;

    // Build rep_info: sample → replicate → "fwd.rev" label
    let mut rep_info: RepInfo = HashMap::new();
    for pool in sample_info.values() {
        for ((ftag, rtag), entry) in pool {
            rep_info
                .entry(entry.sample.clone())
                .or_default()
                .insert(entry.replicate, format!("{}.{}", ftag, rtag));
        }
    }

    // Determine paired-end from first tagInfo header
    let mut paired_end: Option<bool> = None;
    let mut haps_info: HapsInfo = HashMap::new();

    for (pool_name, pool) in &sample_info {
        let path = format!("{}_{}.tagInfo", input_prefix, pool_name);
        let content = fs::read_to_string(&path)
            .with_context(|| format!("Cannot open tagInfo file: {}", path))?;
        let mut lines = content.lines();

        let header_cols = lines.next().unwrap_or("").split('\t').count();
        if paired_end.is_none() {
            paired_end = Some(header_cols == 6);
        }
        let is_paired = paired_end.unwrap();
        let type_index = if is_paired { 5 } else { 4 };

        for line in lines {
            let toks: Vec<&str> = line.split('\t').collect();
            if toks.len() <= type_index || toks[type_index] != "C" {
                continue;
            }
            let tag_pair = (toks[0].to_string(), toks[1].to_string());
            let entry = match pool.get(&tag_pair) {
                Some(e) => e,
                None => continue,
            };
            let (sample, replicate) = (&entry.sample, entry.replicate);
            let (seq, count) = if is_paired {
                (
                    format!("{}_{}", toks[2], toks[3]),
                    toks[4].parse::<u32>().unwrap_or(0),
                )
            } else {
                (toks[2].to_string(), toks[3].parse::<u32>().unwrap_or(0))
            };

            let sample_haps = haps_info.entry(sample.clone()).or_default();
            let seq_entry = sample_haps.entry(seq).or_default();

            // Pre-fill all replicates with 0 on first encounter
            if seq_entry.is_empty() {
                if let Some(reps) = rep_info.get(sample) {
                    for (&rep, label) in reps {
                        seq_entry.insert(rep, (label.clone(), 0));
                    }
                }
            }
            if let Some(slot) = seq_entry.get_mut(&replicate) {
                slot.1 = count;
            }
        }
    }

    let is_paired = paired_end.unwrap_or(false);
    let out_path = format!("{}/{}.fna", output_directory.trim_end_matches('/'), output_prefix);
    let mut out = fs::File::create(&out_path)
        .with_context(|| format!("Cannot create {}", out_path))?;

    for (sample, seqs) in &haps_info {
        for (seq, reps) in seqs {
            let nreps = reps.len();
            let counts: Vec<u32> = reps.values().map(|(_, c)| *c).collect();
            let tag_pairs: Vec<&str> = reps.values().map(|(tp, _)| tp.as_str()).collect();

            let valid = counts.iter().filter(|&&c| c as usize >= min_occurrence).count();
            if (valid as f64 / nreps as f64) < prop_pcrs {
                continue;
            }

            let header = format!(
                ">{}\t{}\t{}",
                sample,
                tag_pairs.join("_"),
                counts.iter().map(|c| c.to_string()).collect::<Vec<_>>().join("_")
            );

            if is_paired {
                let parts: Vec<&str> = seq.splitn(2, '_').collect();
                let r1 = parts.first().copied().unwrap_or("");
                let r2 = parts.get(1).copied().unwrap_or("");
                let combined_len = r1.len() + r2.len();
                if combined_len >= min_length {
                    writeln!(out, "{}\tread1\n{}", header, r1)?;
                    writeln!(out, "{}\tread2\n{}", header, r2)?;
                }
            } else if seq.len() >= min_length {
                writeln!(out, "{}\n{}", header, seq)?;
            }
        }
    }
    Ok(())
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::{tempdir, NamedTempFile};

    fn write_str(path: &std::path::Path, content: &str) {
        std::fs::write(path, content).unwrap();
    }

    #[test]
    fn test_filter_removes_rare_amplicon() {
        // Sample1 has 2 replicates (Pool1 + Pool2).
        // amp_shared (AAATTTGGGCCC) appears in both — should PASS propPCRs=1.0
        // amp_rare (GGGCCCAAATTT) appears only in Pool1 — should be FILTERED
        let dir = tempdir().unwrap();
        let d = dir.path();

        write_str(&d.join("test_Pool1.tagInfo"),
            "FTag\tRTag\tSeq\tCount\tType\ntag1\ttag2\tAAATTTGGGCCC\t2\tC\ntag1\ttag2\tGGGCCCAAATTT\t1\tC\n");
        write_str(&d.join("test_Pool2.tagInfo"),
            "FTag\tRTag\tSeq\tCount\tType\ntag3\ttag4\tAAATTTGGGCCC\t1\tC\n");

        let mut sf = NamedTempFile::new().unwrap();
        writeln!(sf, "Sample1 tag1 tag2 Pool1").unwrap();
        writeln!(sf, "Sample1 tag3 tag4 Pool2").unwrap();

        run_filter(
            d.join("test").to_str().unwrap(),
            sf.path().to_str().unwrap(),
            1.0, 1, 5,
            d.to_str().unwrap(),
            "filtered",
        ).unwrap();

        let fna = fs::read_to_string(d.join("filtered.fna")).unwrap();
        let seqs: Vec<&str> = fna.lines().filter(|l| !l.starts_with('>')).collect();
        assert_eq!(seqs.len(), 1, "only shared amplicon survives propPCRs=1.0");
        assert!(seqs.contains(&"AAATTTGGGCCC"), "shared amplicon must be present");
        assert!(!fna.contains("GGGCCCAAATTT"), "rare amplicon must be filtered");
    }

    #[test]
    fn test_filter_min_occurrence() {
        // Single replicate. amp has count=1 but min_occurrence=2 → filtered.
        let dir = tempdir().unwrap();
        let d = dir.path();

        write_str(&d.join("test_Pool1.tagInfo"),
            "FTag\tRTag\tSeq\tCount\tType\ntag1\ttag2\tAAATTTGGGCCC\t1\tC\n");

        let mut sf = NamedTempFile::new().unwrap();
        writeln!(sf, "Sample1 tag1 tag2 Pool1").unwrap();

        run_filter(
            d.join("test").to_str().unwrap(),
            sf.path().to_str().unwrap(),
            1.0, 2, 5,   // min_occurrence=2
            d.to_str().unwrap(),
            "filtered",
        ).unwrap();

        let fna = fs::read_to_string(d.join("filtered.fna")).unwrap();
        assert!(!fna.contains("AAATTTGGGCCC"), "count=1 < min_occurrence=2 → filtered");
    }

    #[test]
    fn test_filter_min_length() {
        // amp is 5 chars but min_length=10 → filtered
        let dir = tempdir().unwrap();
        let d = dir.path();

        write_str(&d.join("test_Pool1.tagInfo"),
            "FTag\tRTag\tSeq\tCount\tType\ntag1\ttag2\tAAAAA\t5\tC\n");

        let mut sf = NamedTempFile::new().unwrap();
        writeln!(sf, "Sample1 tag1 tag2 Pool1").unwrap();

        run_filter(
            d.join("test").to_str().unwrap(),
            sf.path().to_str().unwrap(),
            1.0, 1, 10,   // min_length=10
            d.to_str().unwrap(),
            "filtered",
        ).unwrap();

        let fna = fs::read_to_string(d.join("filtered.fna")).unwrap();
        assert!(!fna.contains("AAAAA"), "length 5 < min_length 10 → filtered");
    }

    #[test]
    fn test_filter_ignores_non_c_rows() {
        // Row with type B should be ignored
        let dir = tempdir().unwrap();
        let d = dir.path();

        write_str(&d.join("test_Pool1.tagInfo"),
            "FTag\tRTag\tSeq\tCount\tType\ntag1\ttag2\tAAATTTGGGCCC\t2\tC\ntag1\ttag2\tGGGCCCAAATTT\t5\tB\n");

        let mut sf = NamedTempFile::new().unwrap();
        writeln!(sf, "Sample1 tag1 tag2 Pool1").unwrap();

        run_filter(
            d.join("test").to_str().unwrap(),
            sf.path().to_str().unwrap(),
            1.0, 1, 5,
            d.to_str().unwrap(),
            "filtered",
        ).unwrap();

        let fna = fs::read_to_string(d.join("filtered.fna")).unwrap();
        let seqs: Vec<&str> = fna.lines().filter(|l| !l.starts_with('>')).collect();
        assert_eq!(seqs.len(), 1);
        assert!(!fna.contains("GGGCCCAAATTT"), "type B row must not appear in output");
    }
}
