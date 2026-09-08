use std::io;
use std::path::PathBuf;
use std::process::ExitStatus;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum EvidenceError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("genome diff error: {0}")]
    Gd(#[from] prokadiff_gd::GdError),
    #[error("bowtie2 failed ({status}): {stderr}")]
    Bowtie2 { status: ExitStatus, stderr: String },
    #[error("bowtie2-build failed ({status}): {stderr}")]
    Bowtie2Build { status: ExitStatus, stderr: String },
    #[error("required binary not found: {0}")]
    MissingBinary(&'static str),
    #[error("invalid FASTA {path}: {msg}")]
    Fasta { path: PathBuf, msg: String },
    #[error(
        "duplicate reference contig ID {id:?}: first defined in {first_path}, duplicated in {duplicate_path}"
    )]
    DuplicateReferenceId {
        id: String,
        first_path: PathBuf,
        duplicate_path: PathBuf,
    },
    #[error("BAM/SAM error: {0}")]
    Alignment(String),
}

pub type Result<T> = std::result::Result<T, EvidenceError>;
