//! Minimal FASTA / GenBank ORIGIN reader (coordinate skeleton). BAM I/O still uses noodles.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::error::{EvidenceError, Result};

#[derive(Clone, Debug)]
pub struct FastaRecord {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Mobile element or repeat region parsed from reference annotations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepeatRegion {
    pub seq_id: String,
    /// 1-based start coordinate.
    pub start: u64,
    /// 1-based end coordinate.
    pub end: u64,
    pub strand: i8,
    pub name: String,
}

/// Read a FASTA, or a simple GenBank file (LOCUS + ORIGIN) as one record.
pub fn read_reference(path: impl AsRef<Path>) -> Result<Vec<FastaRecord>> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut first = String::new();
    reader.read_line(&mut first)?;
    let trimmed = first.trim_start();
    if trimmed.starts_with('>') {
        drop(reader);
        read_fasta(path)
    } else if trimmed.to_ascii_uppercase().starts_with("LOCUS") {
        drop(reader);
        read_genbank_origin(path)
    } else if trimmed.is_empty() {
        read_fasta(path)
    } else {
        Err(EvidenceError::Fasta {
            path: path.to_path_buf(),
            msg: "expected FASTA (>) or GenBank (LOCUS)".into(),
        })
    }
}

pub fn read_fasta(path: impl AsRef<Path>) -> Result<Vec<FastaRecord>> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut recs = Vec::new();
    let mut name = String::new();
    let mut seq = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if !name.is_empty() {
                recs.push(FastaRecord {
                    name: std::mem::take(&mut name),
                    seq: std::mem::take(&mut seq),
                });
            }
            name = rest.split_whitespace().next().unwrap_or(rest).to_string();
        } else {
            seq.extend(
                line.bytes()
                    .filter(|b| !b.is_ascii_whitespace())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    if !name.is_empty() {
        recs.push(FastaRecord { name, seq });
    }
    if recs.is_empty() {
        return Err(EvidenceError::Fasta {
            path: path.to_path_buf(),
            msg: "no records".into(),
        });
    }
    Ok(recs)
}

fn read_genbank_origin(path: &Path) -> Result<Vec<FastaRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut name = String::from("chr");
    let mut seq = Vec::new();
    let mut in_origin = false;
    for line in reader.lines() {
        let line = line?;
        if !in_origin {
            let t = line.trim_start();
            if t.to_ascii_uppercase().starts_with("LOCUS") {
                if let Some(id) = t.split_whitespace().nth(1) {
                    name = id.to_string();
                }
            } else if t.to_ascii_uppercase().starts_with("ORIGIN") {
                in_origin = true;
            }
            continue;
        }
        if line.trim() == "//" {
            break;
        }
        seq.extend(
            line.bytes()
                .filter(|b| b.is_ascii_alphabetic())
                .map(|b| b.to_ascii_uppercase()),
        );
    }
    if seq.is_empty() {
        return Err(EvidenceError::Fasta {
            path: path.to_path_buf(),
            msg: "GenBank ORIGIN empty".into(),
        });
    }
    Ok(vec![FastaRecord { name, seq }])
}

/// Parse repeat_region and mobile_element annotations from a GenBank file.
pub fn parse_genbank_repeats(path: impl AsRef<Path>) -> Result<Vec<RepeatRegion>> {
    let path = path.as_ref();
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => return Ok(Vec::new()),
    };
    let reader = BufReader::new(file);
    let mut repeats = Vec::new();
    let mut curr_seq = String::from("chr");
    let mut pending_repeat: Option<(u64, u64, i8, String)> = None;

    for line in reader.lines() {
        let line = line?;
        let t = line.trim_start();
        if t.to_ascii_uppercase().starts_with("LOCUS") {
            if let Some(id) = t.split_whitespace().nth(1) {
                curr_seq = id.to_string();
            }
        } else if t.to_ascii_uppercase().starts_with("ORIGIN") {
            if let Some((start, end, strand, name)) = pending_repeat.take() {
                repeats.push(RepeatRegion {
                    seq_id: curr_seq.clone(),
                    start,
                    end,
                    strand,
                    name,
                });
            }
            break;
        } else if line.starts_with("     ") && !line.starts_with("      ") {
            // Feature line (5 spaces indent in standard GenBank)
            if let Some((start, end, strand, name)) = pending_repeat.take() {
                repeats.push(RepeatRegion {
                    seq_id: curr_seq.clone(),
                    start,
                    end,
                    strand,
                    name,
                });
            }
            if t.starts_with("repeat_region") || t.starts_with("mobile_element") {
                let loc_part = t
                    .trim_start_matches("repeat_region")
                    .trim_start_matches("mobile_element")
                    .trim();
                let is_comp = loc_part.starts_with("complement(");
                let s = loc_part
                    .trim_start_matches("complement(")
                    .trim_end_matches(')');
                if let Some((start_s, end_s)) = s.split_once("..") {
                    if let (Ok(start), Ok(end)) =
                        (start_s.trim().parse::<u64>(), end_s.trim().parse::<u64>())
                    {
                        let strand = if is_comp { -1 } else { 1 };
                        pending_repeat = Some((start, end, strand, String::from("repeat")));
                    }
                }
            }
        } else if let Some((_, _, _, ref mut name)) = pending_repeat {
            // Qualifier line for the pending repeat
            let sub_t = line.trim();
            if let Some(rest) = sub_t.strip_prefix("/mobile_element=") {
                let cleaned = rest
                    .trim_matches('"')
                    .trim_start_matches("insertion sequence:")
                    .trim();
                *name = cleaned.to_string();
            } else if *name == "repeat" {
                if let Some(rest) = sub_t.strip_prefix("/note=") {
                    let cleaned = rest.trim_matches('"').trim();
                    *name = cleaned.to_string();
                } else if let Some(rest) = sub_t.strip_prefix("/gene=") {
                    let cleaned = rest.trim_matches('"').trim();
                    *name = cleaned.to_string();
                }
            }
        }
    }
    if let Some((start, end, strand, name)) = pending_repeat.take() {
        repeats.push(RepeatRegion {
            seq_id: curr_seq,
            start,
            end,
            strand,
            name,
        });
    }
    Ok(repeats)
}

/// Concatenate one or more FASTA/GBK paths into a single FASTA file for Bowtie2.
pub fn write_combined_fasta(paths: &[impl AsRef<Path>], dest: &Path) -> Result<()> {
    let mut recs = Vec::new();
    for p in paths {
        recs.extend(read_reference(p.as_ref())?);
    }
    if recs.is_empty() {
        return Err(EvidenceError::Fasta {
            path: dest.to_path_buf(),
            msg: "no reference records".into(),
        });
    }
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)?;
    }
    write_records(&recs, dest)
}

/// Write FASTA records (Bowtie2 index input).
pub fn write_records(recs: &[FastaRecord], dest: &Path) -> Result<()> {
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut out = File::create(dest)?;
    for rec in recs {
        writeln!(out, ">{}", rec.name)?;
        for chunk in rec.seq.chunks(80) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn reads_simple_fasta() {
        let dir = std::env::temp_dir().join("prokdiff-fasta-test");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("t.fa");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, ">chr comment\nACGT\nacgt").unwrap();
        }
        let recs = read_reference(&path).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name, "chr");
        assert_eq!(recs[0].seq, b"ACGTACGT");
    }

    #[test]
    fn reads_genbank_origin() {
        let dir = std::env::temp_dir().join("prokdiff-fasta-test");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("t.gbk");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, "LOCUS       syn  8 bp\nORIGIN\n        1 acgtacgt\n//").unwrap();
        }
        let recs = read_reference(&path).unwrap();
        assert_eq!(recs[0].name, "syn");
        assert_eq!(recs[0].seq, b"ACGTACGT");
    }

    #[test]
    fn parses_genbank_repeat_features() {
        let dir = std::env::temp_dir().join("prokdiff-fasta-test");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("rep.gbk");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(
                f,
                "LOCUS       chr1  1000 bp\nFEATURES             Location/Qualifiers\n     repeat_region   100..200\n                     /mobile_element=\"insertion sequence:IS150\"\n     repeat_region   complement(500..600)\n                     /mobile_element=\"IS186\"\nORIGIN\n        1 aaaa\n//"
            )
            .unwrap();
        }
        let reps = parse_genbank_repeats(&path).unwrap();
        assert_eq!(reps.len(), 2);
        assert_eq!(reps[0].seq_id, "chr1");
        assert_eq!(reps[0].start, 100);
        assert_eq!(reps[0].end, 200);
        assert_eq!(reps[0].strand, 1);
        assert_eq!(reps[0].name, "IS150");

        assert_eq!(reps[1].seq_id, "chr1");
        assert_eq!(reps[1].start, 500);
        assert_eq!(reps[1].end, 600);
        assert_eq!(reps[1].strand, -1);
        assert_eq!(reps[1].name, "IS186");
    }

    #[test]
    fn parses_rel606_gbk_repeats() {
        let path = Path::new("../../testdata/layer2/clonal/Clonal_Sample/REL606.gbk");
        if !path.is_file() {
            return;
        }
        let reps = parse_genbank_repeats(path).unwrap();
        assert!(!reps.is_empty());
        let is150 = reps.iter().filter(|r| r.name == "IS150").count();
        let is186 = reps.iter().filter(|r| r.name == "IS186").count();
        assert!(is150 >= 5, "expected >=5 IS150 copies, found {}", is150);
        assert!(is186 >= 5, "expected >=5 IS186 copies, found {}", is186);
    }
}
