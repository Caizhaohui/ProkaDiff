use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::{EvidenceError, Result};

use super::FastaRecord;

pub(super) fn read_genbank_origin(path: &Path) -> Result<Vec<FastaRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut name = String::from("chr");
    let mut seq = Vec::new();
    let mut in_origin = false;

    for line in reader.lines() {
        let line = line?;
        if !in_origin {
            let trimmed = line.trim_start();
            if trimmed.to_ascii_uppercase().starts_with("LOCUS") {
                if let Some(id) = trimmed.split_whitespace().nth(1) {
                    name = id.to_string();
                }
            } else if trimmed.to_ascii_uppercase().starts_with("ORIGIN") {
                in_origin = true;
            }
            continue;
        }
        if line.trim() == "//" {
            push_record(&mut records, &mut name, &mut seq, path)?;
            in_origin = false;
            continue;
        }
        seq.extend(
            line.bytes()
                .filter(|byte| byte.is_ascii_alphabetic())
                .map(|byte| byte.to_ascii_uppercase()),
        );
    }
    if in_origin {
        push_record(&mut records, &mut name, &mut seq, path)?;
    }
    if records.is_empty() {
        return Err(EvidenceError::Fasta {
            path: path.to_path_buf(),
            msg: "GenBank ORIGIN empty".into(),
        });
    }
    Ok(records)
}

fn push_record(
    records: &mut Vec<FastaRecord>,
    name: &mut String,
    seq: &mut Vec<u8>,
    path: &Path,
) -> Result<()> {
    if seq.is_empty() {
        return Err(EvidenceError::Fasta {
            path: path.to_path_buf(),
            msg: "GenBank ORIGIN empty".into(),
        });
    }
    records.push(FastaRecord {
        name: std::mem::take(name),
        seq: std::mem::take(seq),
    });
    Ok(())
}
