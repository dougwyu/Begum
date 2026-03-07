use crate::SortArgs;
use anyhow::{anyhow, Context, Result};
use indexmap::IndexMap;
use std::{collections::HashMap, fs};

// ── Types ──────────────────────────────────────────────────────────────────

pub struct PoolInfo {
    pub read1: String,
    pub read2: Option<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SampleEntry {
    pub sample: String,
    pub replicate: usize,
}

/// pool_name → { (fwd_tag_name, rev_tag_name) → SampleEntry }
pub type SampleInfo = IndexMap<String, IndexMap<(String, String), SampleEntry>>;

/// tag_name → tag sequence bytes (uppercase)
pub type TagDict = IndexMap<String, Vec<u8>>;

/// (fwd_tag_name, rev_tag_name) → { amplicon_string → count }
pub type Haps = HashMap<(String, String), HashMap<String, u32>>;

// ── File readers ───────────────────────────────────────────────────────────

pub fn read_tag_file(path: &str) -> Result<TagDict> {
    let mut tags = TagDict::new();
    for line in fs::read_to_string(path)
        .with_context(|| format!("Cannot open tag file: {}", path))?
        .lines()
    {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 2 {
            return Err(anyhow!("Tag file line does not have 2 tokens: {}", line));
        }
        tags.insert(toks[0].to_string(), toks[1].to_ascii_uppercase().into_bytes());
    }
    Ok(tags)
}

pub fn read_primer_file(path: &str) -> Result<(Vec<u8>, Vec<u8>)> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Cannot open primer file: {}", path))?;
    let toks: Vec<&str> = content.split_whitespace().collect();
    if toks.len() != 2 {
        return Err(anyhow!(
            "Primer file must have exactly 2 tokens, found {}",
            toks.len()
        ));
    }
    Ok((
        toks[0].to_ascii_uppercase().into_bytes(),
        toks[1].to_ascii_uppercase().into_bytes(),
    ))
}

pub fn read_pool_file(path: &str) -> Result<IndexMap<String, PoolInfo>> {
    let mut pools = IndexMap::new();
    for line in fs::read_to_string(path)
        .with_context(|| format!("Cannot open pool file: {}", path))?
        .lines()
    {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 2 || toks.len() > 3 {
            return Err(anyhow!(
                "Pool file line has wrong number of columns: {}",
                line
            ));
        }
        if pools.contains_key(toks[0]) {
            return Err(anyhow!("Duplicate pool name: {}", toks[0]));
        }
        pools.insert(
            toks[0].to_string(),
            PoolInfo {
                read1: toks[1].to_string(),
                read2: toks.get(2).map(|s| s.to_string()),
            },
        );
    }
    Ok(pools)
}

pub fn read_sample_info_file(path: &str) -> Result<SampleInfo> {
    let mut info: SampleInfo = IndexMap::new();
    let mut sample_count: HashMap<String, usize> = HashMap::new();
    for line in fs::read_to_string(path)
        .with_context(|| format!("Cannot open sample info file: {}", path))?
        .lines()
    {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 4 {
            return Err(anyhow!(
                "Sample info line does not have 4 tokens: {}",
                line
            ));
        }
        let (sample, fwd_tag, rev_tag, pool) = (toks[0], toks[1], toks[2], toks[3]);
        let rep = sample_count.entry(sample.to_string()).or_insert(0);
        *rep += 1;
        let replicate = *rep;
        info.entry(pool.to_string()).or_default().insert(
            (fwd_tag.to_string(), rev_tag.to_string()),
            SampleEntry {
                sample: sample.to_string(),
                replicate,
            },
        );
    }
    Ok(info)
}

/// Find primer positions in two reads.
/// For single-end: pass (seq, RC(seq)). For paired-end: pass (read1, read2).
/// Returns (fstart, fend, rstart, rend, match_type).
/// match_type: 0=none, 1=F-R, 2=F-noR, 3=noF-R, 4=R-F, 5=R-noF, 6=noR-F, 8=F-R-multi, 9=R-F-multi
pub fn find_primer_pos(
    read1: &[u8],
    read2: &[u8],
    fwd_primer: &[u8],
    rev_primer: &[u8],
    max_mismatches: usize,
) -> (i64, i64, i64, i64, u8) {
    use crate::dna;

    // Pass 1: F in read1 (first match), R in read2 (last match)
    let (fs, fe, nf) = dna::find_primer_first(fwd_primer, read1, max_mismatches);
    let (rs, re, nr) = dna::find_primer_last(rev_primer, read2, max_mismatches);

    if fs != -1 || rs != -1 {
        let mt = if fs != -1 && rs != -1 {
            if nf > 1 || nr > 1 { 8 } else { 1 }
        } else if fs != -1 {
            2
        } else {
            3
        };
        return (fs, fe, rs, re, mt);
    }

    // Pass 2: R in read1 (first match), F in read2 (last match)
    let (rs2, re2, nr2) = dna::find_primer_first(rev_primer, read1, max_mismatches);
    let (fs2, fe2, nf2) = dna::find_primer_last(fwd_primer, read2, max_mismatches);

    let mt = if fs2 != -1 && rs2 != -1 {
        if nf2 > 1 || nr2 > 1 { 9 } else { 4 }
    } else if rs2 != -1 {
        5
    } else if fs2 != -1 {
        6
    } else {
        0
    };
    (fs2, fe2, rs2, re2, mt)
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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

pub fn find_best_tag_match(
    tag_region: &[u8],
    tags: &TagDict,
    max_errors: usize,
) -> Option<String> {
    let mut best_name: Option<String> = None;
    let mut best_dist = usize::MAX;
    let mut num_matches = 0usize;

    for (name, seq) in tags {
        let dist = crate::dna::hamming_tag_distance(seq, tag_region);
        if dist < best_dist {
            best_dist = dist;
            best_name = Some(name.clone());
            num_matches = 1;
        } else if dist == best_dist {
            // Tie-break: longer tag wins (mirrors Python)
            let best_len = best_name.as_ref().map(|n| tags[n].len()).unwrap_or(0);
            if seq.len() > best_len {
                best_name = Some(name.clone());
                num_matches = 1;
            } else if seq.len() == best_len {
                num_matches += 1;
            }
        }
    }

    if best_dist > max_errors || num_matches > 1 {
        return None;
    }
    best_name
}

pub fn process_single_end(
    path: &str,
    fwd_primer: &[u8],
    rev_primer: &[u8],
    tags: &TagDict,
    primer_mismatches: usize,
    tag_mismatches: usize,
    allow_multi_primers: bool,
    primer_counts: &mut [i64; 10],
    tag_counts: &mut [i64; 8],
) -> Result<Haps> {
    let mut haps: Haps = HashMap::new();
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("Cannot open FASTQ: {}", path))?;

    while let Some(rec) = reader.next() {
        let rec = rec?;
        let seq = rec.seq();
        let seq_rc = reverse_complement(&seq);

        let (fs, fe, rs, re, mt) =
            find_primer_pos(&seq, &seq_rc, fwd_primer, rev_primer, primer_mismatches);
        primer_counts[mt as usize] += 1;

        let (ftag_str, rtag_str, amplicon) =
            if (allow_multi_primers && mt == 8) || mt == 1 {
                let fs = fs as usize;
                let fe = fe as usize;
                let rs = rs as usize;
                let re = re as usize;
                let seq_len = seq.len();
                (
                    String::from_utf8_lossy(&seq[..fs]).to_string(),
                    String::from_utf8_lossy(&seq_rc[..rs]).to_string(),
                    String::from_utf8_lossy(&seq[fe..seq_len - re]).to_string(),
                )
            } else if (allow_multi_primers && mt == 9) || mt == 4 {
                let fs = fs as usize;
                let fe = fe as usize;
                let rs = rs as usize;
                let re = re as usize;
                let seq_len = seq.len();
                (
                    String::from_utf8_lossy(&seq_rc[..fs]).to_string(),
                    String::from_utf8_lossy(&seq[..rs]).to_string(),
                    String::from_utf8_lossy(&seq_rc[fe..seq_len - re]).to_string(),
                )
            } else {
                continue;
            };

        if amplicon.is_empty() {
            primer_counts[7] += 1;
            continue;
        }

        let best_ftag = find_best_tag_match(ftag_str.as_bytes(), tags, tag_mismatches);
        let best_rtag = find_best_tag_match(rtag_str.as_bytes(), tags, tag_mismatches);
        let tag_type =
            (best_ftag.is_none()) as usize + 2 * (best_rtag.is_none()) as usize;
        let orient_offset = if mt == 4 { 4 } else { 0 };
        tag_counts[orient_offset + tag_type] += 1;

        if tag_type != 0 {
            continue;
        }

        let tag_combo = (best_ftag.unwrap(), best_rtag.unwrap());
        *haps
            .entry(tag_combo)
            .or_default()
            .entry(amplicon)
            .or_insert(0) += 1;
    }
    Ok(haps)
}

pub fn process_paired_end(
    path1: &str,
    path2: &str,
    fwd_primer: &[u8],
    rev_primer: &[u8],
    tags: &TagDict,
    primer_mismatches: usize,
    tag_mismatches: usize,
    primer_counts: &mut [i64; 10],
    tag_counts: &mut [i64; 8],
) -> Result<Haps> {
    let mut haps: Haps = HashMap::new();
    let mut r1 = needletail::parse_fastx_file(path1)
        .with_context(|| format!("Cannot open {}", path1))?;
    let mut r2 = needletail::parse_fastx_file(path2)
        .with_context(|| format!("Cannot open {}", path2))?;

    loop {
        match (r1.next(), r2.next()) {
            (Some(rec1), Some(rec2)) => {
                let rec1 = rec1?;
                let rec2 = rec2?;
                let seq1 = rec1.seq();
                let seq2 = rec2.seq();

                let (fs, fe, rs, re, mt) = find_primer_pos(
                    &seq1, &seq2, fwd_primer, rev_primer, primer_mismatches,
                );
                primer_counts[mt as usize] += 1;

                let (ftag_str, rtag_str, amplicon) = if mt == 1 {
                    (
                        String::from_utf8_lossy(&seq1[..fs as usize]).to_string(),
                        String::from_utf8_lossy(&seq2[..rs as usize]).to_string(),
                        format!(
                            "{}\t{}",
                            String::from_utf8_lossy(&seq1[fe as usize..]),
                            String::from_utf8_lossy(&seq2[re as usize..])
                        ),
                    )
                } else if mt == 4 {
                    (
                        String::from_utf8_lossy(&seq2[..rs as usize]).to_string(),
                        String::from_utf8_lossy(&seq1[..fs as usize]).to_string(),
                        format!(
                            "{}\t{}",
                            String::from_utf8_lossy(&seq2[re as usize..]),
                            String::from_utf8_lossy(&seq1[fe as usize..])
                        ),
                    )
                } else {
                    continue;
                };

                let best_ftag =
                    find_best_tag_match(ftag_str.as_bytes(), tags, tag_mismatches);
                let best_rtag =
                    find_best_tag_match(rtag_str.as_bytes(), tags, tag_mismatches);
                let tag_type =
                    (best_ftag.is_none()) as usize + 2 * (best_rtag.is_none()) as usize;
                let orient_offset = if mt == 4 { 4 } else { 0 };
                tag_counts[orient_offset + tag_type] += 1;

                if tag_type != 0 {
                    continue;
                }

                let tag_combo = (best_ftag.unwrap(), best_rtag.unwrap());
                *haps
                    .entry(tag_combo)
                    .or_default()
                    .entry(amplicon)
                    .or_insert(0) += 1;
            }
            (None, None) => break,
            _ => return Err(anyhow!("Paired FASTQ files have different read counts")),
        }
    }
    Ok(haps)
}

pub fn write_tag_info(
    haps: &Haps,
    outprefix: &str,
    pool_name: &str,
    sample_info: &SampleInfo,
    single_end: bool,
) -> Result<()> {
    use std::io::Write;
    let path = format!("{}_{}.tagInfo", outprefix, pool_name);
    let mut f = fs::File::create(&path)
        .with_context(|| format!("Cannot create {}", path))?;

    if single_end {
        writeln!(f, "FTag\tRTag\tSeq\tCount\tType")?;
    } else {
        writeln!(f, "FTag\tRTag\tFSeq\tRSeq\tCount\tType")?;
    }

    // Collect pool tags and valid pairs
    let pool_tags: std::collections::HashSet<String> = sample_info
        .get(pool_name)
        .map(|p| {
            p.keys()
                .flat_map(|(a, b)| [a.clone(), b.clone()])
                .collect()
        })
        .unwrap_or_default();
    let pool_map = sample_info.get(pool_name);

    for ((ftag, rtag), seqs) in haps {
        let ftag_in = pool_tags.contains(ftag);
        let rtag_in = pool_tags.contains(rtag);
        let tag_type = if ftag_in && rtag_in {
            if pool_map.map_or(false, |p| p.contains_key(&(ftag.clone(), rtag.clone()))) {
                "C"
            } else {
                "B"
            }
        } else if ftag_in {
            "F"
        } else if rtag_in {
            "R"
        } else {
            "N"
        };

        for (seq, count) in seqs {
            if single_end {
                writeln!(f, "{}\t{}\t{}\t{}\t{}", ftag, rtag, seq, count, tag_type)?;
            } else {
                let parts: Vec<&str> = seq.splitn(2, '\t').collect();
                writeln!(
                    f,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    ftag,
                    rtag,
                    parts.first().copied().unwrap_or(""),
                    parts.get(1).copied().unwrap_or(""),
                    count,
                    tag_type
                )?;
            }
        }
    }
    Ok(())
}

pub fn run(args: SortArgs) -> Result<()> {
    let tags = read_tag_file(&args.tags)?;

    let (fwd_primer, rev_primer) = if let Some(pf) = &args.primers {
        read_primer_file(pf)?
    } else {
        match (&args.fwd_primer, &args.rev_primer) {
            (Some(f), Some(r)) => (
                f.to_ascii_uppercase().into_bytes(),
                r.to_ascii_uppercase().into_bytes(),
            ),
            _ => {
                return Err(anyhow!(
                    "Specify a primer file (-p) or both --p1 and --p2"
                ))
            }
        }
    };

    let pools = read_pool_file(&args.pool)?;
    let sample_info = read_sample_info_file(&args.sample_info)?;

    let out_dir = args.output_directory.trim_end_matches('/');
    let outprefix = if args.output_prefix.is_empty() {
        out_dir.to_string()
    } else {
        format!("{}/{}", out_dir, args.output_prefix)
    };

    for (pool_name, pool) in &pools {
        log::info!("Processing pool {}", pool_name);
        let mut primer_counts = [0i64; 10];
        let mut tag_counts = [0i64; 8];

        let haps = if pool.read2.is_none() {
            log::info!("Single-end mode for pool {}", pool_name);
            process_single_end(
                &pool.read1,
                &fwd_primer,
                &rev_primer,
                &tags,
                args.primer_mismatches,
                args.tag_mismatches,
                args.allow_multiple_primers,
                &mut primer_counts,
                &mut tag_counts,
            )?
        } else {
            log::info!("Paired-end mode for pool {}", pool_name);
            process_paired_end(
                &pool.read1,
                pool.read2.as_deref().unwrap(),
                &fwd_primer,
                &rev_primer,
                &tags,
                args.primer_mismatches,
                args.tag_mismatches,
                &mut primer_counts,
                &mut tag_counts,
            )?
        };

        let single_end = pool.read2.is_none();
        write_tag_info(&haps, &outprefix, pool_name, &sample_info, single_end)?;

        log::info!(
            "Pool {} done. Primer match counts: {:?}",
            pool_name,
            primer_counts
        );
        log::info!("Tag match counts: {:?}", tag_counts);
    }

    log::info!("Sorting complete.");
    Ok(())
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn test_read_tag_file_basic() {
        let f = write_temp("tag1 AACCGGT\ntag2 TTGGCCA\n");
        let tags = read_tag_file(f.path().to_str().unwrap()).unwrap();
        assert_eq!(tags.len(), 2);
        assert_eq!(tags["tag1"], b"AACCGGT".to_vec());
        assert_eq!(tags["tag2"], b"TTGGCCA".to_vec());
    }

    #[test]
    fn test_read_tag_file_uppercase() {
        // sequences must be stored uppercase regardless of input
        let f = write_temp("tag1 aaccggt\n");
        let tags = read_tag_file(f.path().to_str().unwrap()).unwrap();
        assert_eq!(tags["tag1"], b"AACCGGT".to_vec());
    }

    #[test]
    fn test_read_tag_file_skips_blank_lines() {
        let f = write_temp("tag1 AACCGGT\n\ntag2 TTGGCCA\n");
        let tags = read_tag_file(f.path().to_str().unwrap()).unwrap();
        assert_eq!(tags.len(), 2);
    }

    #[test]
    fn test_read_primer_file() {
        let f = write_temp("GCATGC AGTCAG\n");
        let (fwd, rev) = read_primer_file(f.path().to_str().unwrap()).unwrap();
        assert_eq!(fwd, b"GCATGC".to_vec());
        assert_eq!(rev, b"AGTCAG".to_vec());
    }

    #[test]
    fn test_read_pool_file_single_end() {
        let f = write_temp("Pool1 reads.fastq\n");
        let pools = read_pool_file(f.path().to_str().unwrap()).unwrap();
        assert_eq!(pools.len(), 1);
        assert_eq!(pools["Pool1"].read1, "reads.fastq");
        assert!(pools["Pool1"].read2.is_none());
    }

    #[test]
    fn test_read_pool_file_paired_end() {
        let f = write_temp("Pool1 r1.fastq r2.fastq\n");
        let pools = read_pool_file(f.path().to_str().unwrap()).unwrap();
        assert_eq!(pools["Pool1"].read2, Some("r2.fastq".to_string()));
    }

    #[test]
    fn test_read_pool_file_multiple_pools() {
        let f = write_temp("Pool1 r1.fastq\nPool2 r2.fastq\n");
        let pools = read_pool_file(f.path().to_str().unwrap()).unwrap();
        assert_eq!(pools.len(), 2);
    }

    #[test]
    fn test_read_pool_file_duplicate_name_errors() {
        let f = write_temp("Pool1 r1.fastq\nPool1 r2.fastq\n");
        assert!(read_pool_file(f.path().to_str().unwrap()).is_err());
    }

    #[test]
    fn test_read_sample_info_replicates() {
        // Sample1 appears twice → replicate 1 and 2
        let f = write_temp("Sample1 tag1 tag2 Pool1\nSample1 tag3 tag4 Pool2\n");
        let info = read_sample_info_file(f.path().to_str().unwrap()).unwrap();

        let entry1 = &info["Pool1"][&("tag1".to_string(), "tag2".to_string())];
        assert_eq!(entry1.sample, "Sample1");
        assert_eq!(entry1.replicate, 1);

        let entry2 = &info["Pool2"][&("tag3".to_string(), "tag4".to_string())];
        assert_eq!(entry2.sample, "Sample1");
        assert_eq!(entry2.replicate, 2);
    }

    #[test]
    fn test_read_sample_info_different_samples() {
        let f = write_temp("SampleA tag1 tag2 Pool1\nSampleB tag3 tag4 Pool1\n");
        let info = read_sample_info_file(f.path().to_str().unwrap()).unwrap();
        let ea = &info["Pool1"][&("tag1".to_string(), "tag2".to_string())];
        let eb = &info["Pool1"][&("tag3".to_string(), "tag4".to_string())];
        assert_eq!(ea.replicate, 1); // SampleA first occurrence
        assert_eq!(eb.replicate, 1); // SampleB first occurrence
    }

    // ── reverse_complement tests ───────────────────────────────────────────

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement(b"GCATGC"), b"GCATGC"); // palindrome
        assert_eq!(reverse_complement(b"AACCGGT"), b"ACCGGTT");
        assert_eq!(reverse_complement(b"AGTCAG"), b"CTGACT");
    }

    // ── find_best_tag_match tests ──────────────────────────────────────────

    #[test]
    fn test_find_best_tag_match_exact() {
        let mut tags = TagDict::new();
        tags.insert("tag1".into(), b"AACCGGT".to_vec());
        tags.insert("tag2".into(), b"TTGGCCA".to_vec());
        assert_eq!(
            find_best_tag_match(b"AACCGGT", &tags, 0),
            Some("tag1".to_string())
        );
    }

    #[test]
    fn test_find_best_tag_match_suffix() {
        // tag region may be longer than tag — compares against end
        let mut tags = TagDict::new();
        tags.insert("tag1".into(), b"AACCGGT".to_vec());
        // "XXXAACCGGT" — last 7 chars match tag1
        assert_eq!(
            find_best_tag_match(b"XXXAACCGGT", &tags, 0),
            Some("tag1".to_string())
        );
    }

    #[test]
    fn test_find_best_tag_match_no_match() {
        let mut tags = TagDict::new();
        tags.insert("tag1".into(), b"AACCGGT".to_vec());
        assert_eq!(find_best_tag_match(b"TTTTTTT", &tags, 0), None);
    }

    #[test]
    fn test_find_best_tag_match_within_errors() {
        let mut tags = TagDict::new();
        tags.insert("tag1".into(), b"AACCGGT".to_vec());
        // 1 mismatch: AACCGGG vs AACCGGT
        assert_eq!(
            find_best_tag_match(b"AACCGGG", &tags, 1),
            Some("tag1".to_string())
        );
        // same mismatch but errors=0 → no match
        assert_eq!(find_best_tag_match(b"AACCGGG", &tags, 0), None);
    }

    // ── process_single_end tests ───────────────────────────────────────────

    #[test]
    fn test_process_single_end_basic() {
        // Construct read: tag1 + fwd_primer + amplicon + RC(rev_primer) + RC(tag2)
        // tag1=AACCGGT, tag2=TTGGCCA, fwd=GCATGC, rev=AGTCAG
        // RC(AGTCAG) = CTGACT, RC(TTGGCCA) = TGGCCAA
        let seq = b"AACCGGTGCATGCAAATTTGGGCCCCTGACTTGGCCAA";
        let qual = "I".repeat(seq.len());
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "@read1\n{}\n+\n{}", std::str::from_utf8(seq).unwrap(), qual).unwrap();
        writeln!(f, "@read2\n{}\n+\n{}", std::str::from_utf8(seq).unwrap(), qual).unwrap();

        let mut tags = TagDict::new();
        tags.insert("tag1".into(), b"AACCGGT".to_vec());
        tags.insert("tag2".into(), b"TTGGCCA".to_vec());

        let mut primer_counts = [0i64; 10];
        let mut tag_counts = [0i64; 8];

        let haps = process_single_end(
            f.path().to_str().unwrap(),
            b"GCATGC", b"AGTCAG",
            &tags, 0, 0, false,
            &mut primer_counts, &mut tag_counts,
        ).unwrap();

        let key = ("tag1".to_string(), "tag2".to_string());
        assert!(haps.contains_key(&key), "tag1/tag2 combination should be found");
        assert_eq!(haps[&key]["AAATTTGGGCCC"], 2, "amplicon seen twice");
        assert_eq!(primer_counts[1], 2, "2 F-R matches");
        assert_eq!(tag_counts[0], 2, "2 both-tags-found");
    }

    #[test]
    fn test_process_single_end_no_primer() {
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let qual = "I".repeat(seq.len());
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "@read1\n{}\n+\n{}", std::str::from_utf8(seq).unwrap(), qual).unwrap();

        let mut tags = TagDict::new();
        tags.insert("tag1".into(), b"AACCGGT".to_vec());

        let mut primer_counts = [0i64; 10];
        let mut tag_counts = [0i64; 8];

        let haps = process_single_end(
            f.path().to_str().unwrap(),
            b"GCATGC", b"AGTCAG",
            &tags, 0, 0, false,
            &mut primer_counts, &mut tag_counts,
        ).unwrap();

        assert!(haps.is_empty());
        assert_eq!(primer_counts[0], 1, "1 no-match read");
    }

    // ── find_primer_pos tests ─────────────────────────────────────────────

    #[test]
    fn test_find_primer_pos_fr_match() {
        // Read: [tag1=AACCGGT][fwd=GCATGC][amplicon][rc_rev=CTGACT][rc_tag2=TGGCCAA]
        // seq    = AACCGGT GCATGC AAATTTGGGCCC CTGACT TGGCCAA
        // seq_rc = TTGGCCAA GTCAG GGGCCCAAATTT GCATGC ACCGGTT
        // (we only need fwd_primer found in seq, rev_primer found in seq_rc)
        let fwd = b"GCATGC";
        let rev = b"AGTCAG";
        let seq:    &[u8] = b"AACCGGTGCATGCAAATTTGGGCCCCTGACTTGGCCAA";
        // rc_rev = CTGACT, rc_tag2 = TGGCCAA
        // For seq_rc we need AGTCAG to appear — it appears in the RC of seq
        // RC of seq: TTGGCCAAGTCAGGGGCCCAAATTTGCATGCACCGGTT
        let seq_rc: &[u8] = b"TTGGCCAAGTCAGGGGCCCAAATTTGCATGCACCGGTT";
        let (fs, fe, rs, re, mt) = find_primer_pos(seq, seq_rc, fwd, rev, 0);
        assert_eq!(mt, 1, "should be F-R match type 1");
        assert_eq!(fs, 7,  "fwd primer starts at pos 7 (after 7-char tag)");
        assert_eq!(fe, 13, "fwd primer ends at pos 13");
        assert!(rs >= 0,   "rev primer found in seq_rc");
        assert!(re > rs,   "rev primer end > start");
    }

    #[test]
    fn test_find_primer_pos_no_match() {
        let fwd = b"GCATGC";
        let rev = b"AGTCAG";
        let seq:    &[u8] = b"AAAAAAAAAAAAAAAAAAA";
        let seq_rc: &[u8] = b"TTTTTTTTTTTTTTTTTTT";
        let (_fs, _fe, _rs, _re, mt) = find_primer_pos(seq, seq_rc, fwd, rev, 0);
        assert_eq!(mt, 0, "should be no match");
    }

    #[test]
    fn test_find_primer_pos_f_only() {
        // F primer found in read1 but R not found in read2
        let fwd = b"GCATGC";
        let rev = b"AGTCAG";
        let read1: &[u8] = b"AACCGGTGCATGCAAATTT";  // has fwd
        let read2: &[u8] = b"TTTTTTTTTTTTTTTTTTT";  // no rev
        let (_fs, _fe, _rs, _re, mt) = find_primer_pos(read1, read2, fwd, rev, 0);
        assert_eq!(mt, 2, "F found, R not found = type 2");
    }

    #[test]
    fn test_find_primer_pos_rf_match() {
        // R primer in read1, F primer in read2 (reverse orientation, type 4)
        let fwd = b"GCATGC";
        let rev = b"AGTCAG";
        let read1: &[u8] = b"AACCGGTAGTCAGAAATTT";  // has rev primer
        // read2 needs fwd at end (find_last_match)
        let read2b: &[u8] = b"TTGGCCAGCATGCZZZ";
        let (_fs, _fe, _rs, _re, mt) = find_primer_pos(read1, read2b, fwd, rev, 0);
        assert_eq!(mt, 4, "R in read1, F in read2 = type 4");
    }

    // ── write_tag_info tests ───────────────────────────────────────────────

    #[test]
    fn test_write_tag_info_single_end() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("test").to_str().unwrap().to_string();

        let mut haps: Haps = HashMap::new();
        haps.entry(("tag1".to_string(), "tag2".to_string()))
            .or_default()
            .insert("AAATTTGGGCCC".to_string(), 2);
        haps.entry(("tag1".to_string(), "tag2".to_string()))
            .or_default()
            .insert("GGGCCCAAATTT".to_string(), 1);
        // unknown tag combo (both not in pool)
        haps.entry(("tagX".to_string(), "tagY".to_string()))
            .or_default()
            .insert("CCCCCCCCCCCC".to_string(), 1);

        let mut sample_info: SampleInfo = IndexMap::new();
        sample_info
            .entry("Pool1".to_string())
            .or_default()
            .insert(
                ("tag1".to_string(), "tag2".to_string()),
                SampleEntry { sample: "Sample1".to_string(), replicate: 1 },
            );

        write_tag_info(&haps, &prefix, "Pool1", &sample_info, true).unwrap();

        let content = std::fs::read_to_string(format!("{}_Pool1.tagInfo", prefix)).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        // Header
        assert_eq!(lines[0], "FTag\tRTag\tSeq\tCount\tType");
        // tag1/tag2 rows should be type C
        let c_rows: Vec<&&str> = lines.iter().filter(|l| l.ends_with("\tC")).collect();
        assert_eq!(c_rows.len(), 2, "2 correct-combo rows");
        // tagX/tagY row should be type N
        let n_rows: Vec<&&str> = lines.iter().filter(|l| l.ends_with("\tN")).collect();
        assert_eq!(n_rows.len(), 1, "1 neither-tag row");
        // Specific count check
        assert!(content.contains("tag1\ttag2\tAAATTTGGGCCC\t2\tC"));
        assert!(content.contains("tag1\ttag2\tGGGCCCAAATTT\t1\tC"));
    }
}
