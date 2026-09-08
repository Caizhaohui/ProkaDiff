//! Minimal FASTA / GenBank ORIGIN reader (coordinate skeleton). BAM I/O still uses noodles.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::error::{EvidenceError, Result};

mod genbank;
mod reference_validation;

use genbank::read_genbank_origin;
use reference_validation::validate_reference_ids;

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

/// Read FASTA or GenBank LOCUS/ORIGIN records from one reference source.
pub fn read_reference(path: impl AsRef<Path>) -> Result<Vec<FastaRecord>> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut first = String::new();
    reader.read_line(&mut first)?;
    let trimmed = first.trim_start();
    let records = if trimmed.starts_with('>') {
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
    }?;
    let mut ids = Default::default();
    validate_reference_ids(&records, path, &mut ids)?;
    Ok(records)
}

/// Read one or more reference inputs after checking that contig IDs are globally unique.
pub fn read_references(paths: &[impl AsRef<Path>]) -> Result<Vec<FastaRecord>> {
    let mut ids = Default::default();
    let mut records = Vec::new();
    for path in paths {
        let source = path.as_ref();
        let source_records = read_reference(source)?;
        validate_reference_ids(&source_records, source, &mut ids)?;
        records.extend(source_records);
    }
    Ok(records)
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
    let recs = read_references(paths)?;
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
mod tests;
