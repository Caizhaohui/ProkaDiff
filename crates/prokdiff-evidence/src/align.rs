//! External Bowtie2 wrapper. Converts SAM → coordinate-sorted BAM with noodles.

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as _;

use crate::error::{EvidenceError, Result};

#[derive(Clone, Debug)]
pub struct FastqInput {
    pub files: Vec<PathBuf>,
}

impl FastqInput {
    pub fn se(path: impl Into<PathBuf>) -> Self {
        Self {
            files: vec![path.into()],
        }
    }

    pub fn pe(r1: impl Into<PathBuf>, r2: impl Into<PathBuf>) -> Self {
        Self {
            files: vec![r1.into(), r2.into()],
        }
    }
}

pub fn which_or_err(name: &'static str) -> Result<PathBuf> {
    which(name).ok_or(EvidenceError::MissingBinary(name))
}

fn which(name: &str) -> Option<PathBuf> {
    if let Ok(p) = std::env::var(format!("{}_PATH", name.to_ascii_uppercase())) {
        let pb = PathBuf::from(p);
        if pb.is_file() {
            return Some(pb);
        }
    }
    std::env::var_os("PATH").and_then(|paths| {
        std::env::split_paths(&paths).find_map(|dir| {
            let p = dir.join(name);
            p.is_file().then_some(p)
        })
    })
}

/// `bowtie2-build` + `bowtie2 --local`, then noodles SAM→sorted BAM.
pub fn align_to_bam(
    ref_fa: &Path,
    reads: &FastqInput,
    out_bam: &Path,
    threads: usize,
    work: &Path,
) -> Result<()> {
    let bowtie2 = which_or_err("bowtie2")?;
    let build = which("bowtie2-build").unwrap_or_else(|| {
        bowtie2
            .parent()
            .map(|d| d.join("bowtie2-build"))
            .unwrap_or_else(|| PathBuf::from("bowtie2-build"))
    });
    if !build.is_file() && which("bowtie2-build").is_none() {
        return Err(EvidenceError::MissingBinary("bowtie2-build"));
    }

    std::fs::create_dir_all(work)?;
    let prefix = work.join("bt2ref");
    let status = Command::new(&build)
        .args([
            "--threads",
            &threads.max(1).to_string(),
            ref_fa.to_str().unwrap_or(""),
            prefix.to_str().unwrap_or(""),
        ])
        .output()?;
    if !status.status.success() {
        return Err(EvidenceError::Bowtie2Build {
            status: status.status,
            stderr: String::from_utf8_lossy(&status.stderr).into_owned(),
        });
    }

    let sam_path = work.join("aligned.sam");
    let mut cmd = Command::new(&bowtie2);
    cmd.args([
        "--local",
        "--no-unal",
        "-p",
        &threads.max(1).to_string(),
        "-x",
        prefix.to_str().unwrap_or(""),
        "-S",
        sam_path.to_str().unwrap_or(""),
    ]);
    match reads.files.as_slice() {
        [se] => {
            cmd.arg("-U").arg(se);
        }
        [r1, r2] => {
            cmd.arg("-1").arg(r1).arg("-2").arg(r2);
        }
        files if files.len() >= 2 && files.len() % 2 == 0 => {
            let mut r1: Vec<String> = Vec::new();
            let mut r2: Vec<String> = Vec::new();
            for chunk in files.chunks(2) {
                r1.push(chunk[0].to_string_lossy().into_owned());
                r2.push(chunk[1].to_string_lossy().into_owned());
            }
            cmd.arg("-1").arg(r1.join(",")).arg("-2").arg(r2.join(","));
        }
        extra => {
            return Err(EvidenceError::Alignment(format!(
                "expected 1 FASTQ (SE) or an even number of FASTQ files (PE pairs), got {}",
                extra.len()
            )));
        }
    }
    let out = cmd.output()?;
    if !out.status.success() {
        return Err(EvidenceError::Bowtie2 {
            status: out.status,
            stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
        });
    }

    sam_to_sorted_bam(&sam_path, out_bam)?;
    Ok(())
}

fn sam_to_sorted_bam(sam_path: &Path, bam_path: &Path) -> Result<()> {
    let mut reader = File::open(sam_path)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let mut records: Vec<sam::Record> = Vec::new();
    for rec in reader.records() {
        records.push(rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?);
    }
    records.sort_by(|a, b| {
        let aid = a
            .reference_sequence_id(&header)
            .and_then(std::result::Result::ok);
        let bid = b
            .reference_sequence_id(&header)
            .and_then(std::result::Result::ok);
        let apos = a
            .alignment_start()
            .and_then(std::result::Result::ok)
            .map(|p| p.get())
            .unwrap_or(0);
        let bpos = b
            .alignment_start()
            .and_then(std::result::Result::ok)
            .map(|p| p.get())
            .unwrap_or(0);
        aid.cmp(&bid).then(apos.cmp(&bpos))
    });

    if let Some(parent) = bam_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut writer = File::create(bam_path).map(bam::io::Writer::new)?;
    writer
        .write_header(&header)
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    for rec in &records {
        writer
            .write_alignment_record(&header, rec)
            .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    }
    writer
        .try_finish()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pe_input_keeps_r1_then_r2_order() {
        let fq = FastqInput::pe("a_R1.fastq", "a_R2.fastq");
        assert_eq!(fq.files.len(), 2);
        assert!(fq.files[0].ends_with("a_R1.fastq"));
        assert!(fq.files[1].ends_with("a_R2.fastq"));
    }
}
