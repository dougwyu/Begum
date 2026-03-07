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

pub fn run(_args: SortArgs) -> Result<()> {
    todo!("sort not yet implemented")
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
        let read2: &[u8] = b"TTGGCCAGCATGCAAATTT"; // has fwd primer (GCATGC at pos 7... wait)
        // read2 needs fwd at end (find_last_match). Let's put fwd at the end.
        // read2 = TTGGCCA + GCATGC + ...  => find_last finds it
        let read2b: &[u8] = b"TTGGCCAGCATGCZZZ";
        let (_fs, _fe, _rs, _re, mt) = find_primer_pos(read1, read2b, fwd, rev, 0);
        assert_eq!(mt, 4, "R in read1, F in read2 = type 4");
    }
}
