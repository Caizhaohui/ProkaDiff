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
}
