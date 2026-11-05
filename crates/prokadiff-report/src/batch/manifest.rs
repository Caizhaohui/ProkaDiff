//! Multi-clone batch manifest parsing and validation.

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use thiserror::Error;

/// An individual clone entry defined in a multi-clone batch manifest.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BatchManifestEntry {
    /// Unique alphanumeric identifier for the clone.
    pub clone_id: String,
    /// Path to R1 FASTQ.
    pub r1: PathBuf,
    /// Optional path to R2 FASTQ for paired-end sequencing.
    pub r2: Option<PathBuf>,
    /// Optional per-clone intended edit table path (overrides global `--intended`).
    pub intended_tsv: Option<PathBuf>,
    /// Optional per-clone spacer sequence override.
    pub spacer: Option<String>,
    /// Optional per-clone PAM profile override.
    pub pam: Option<String>,
}

/// A parsed and validated multi-clone batch manifest.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BatchManifest {
    pub entries: Vec<BatchManifestEntry>,
}

#[derive(Debug, Error, PartialEq, Eq)]
pub enum BatchManifestError {
    #[error("batch manifest file not found: {}", .0.display())]
    FileNotFound(PathBuf),
    #[error("failed to read batch manifest: {0}")]
    Io(String),
    #[error("batch manifest is empty")]
    EmptyManifest,
    #[error("batch manifest missing required header columns: {0}")]
    MissingRequiredColumns(String),
    #[error("line {line}: ragged row (expected {expected} columns, got {actual})")]
    RaggedRow {
        line: usize,
        expected: usize,
        actual: usize,
    },
    #[error("line {line}: clone_id is empty")]
    EmptyCloneId { line: usize },
    #[error("line {line}: invalid clone_id '{clone_id}': {reason}")]
    InvalidCloneId {
        line: usize,
        clone_id: String,
        reason: String,
    },
    #[error("line {line}: duplicate clone_id '{clone_id}'")]
    DuplicateCloneId { line: usize, clone_id: String },
    #[error("line {line}: r1 path is empty for clone '{clone_id}'")]
    EmptyR1 { line: usize, clone_id: String },
    #[error("clone '{clone_id}': {field} file not found: {}", path.display())]
    InputFileNotFound {
        clone_id: String,
        path: PathBuf,
        field: &'static str,
    },
    #[error("line {line}: clone '{clone_id}': invalid spacer sequence: {reason}")]
    InvalidSpacer {
        line: usize,
        clone_id: String,
        reason: String,
    },
    #[error("line {line}: clone '{clone_id}': invalid PAM motif: {reason}")]
    InvalidPam {
        line: usize,
        clone_id: String,
        reason: String,
    },
}

/// Validate that a clone ID conforms to safe, single-component identifier grammar.
/// Must be alphanumeric, `_`, or `-`, length between 1 and 64, no path separators or directory traversals.
pub fn validate_clone_id(raw: &str) -> Result<String, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("clone_id cannot be empty".to_string());
    }
    if trimmed.len() > 64 {
        return Err("clone_id exceeds maximum length of 64 characters".to_string());
    }
    if trimmed == "." || trimmed == ".." {
        return Err("clone_id cannot be '.' or '..'".to_string());
    }
    for c in trimmed.chars() {
        if !c.is_ascii_alphanumeric() && c != '_' && c != '-' {
            return Err(format!(
                "contains invalid character '{c}'; expected alphanumeric, '_' or '-'"
            ));
        }
    }
    Ok(trimmed.to_string())
}

/// Validate spacer sequence: non-empty, DNA alphabet (no U/u), A/C/G/T/N only.
pub fn validate_manifest_spacer(raw: &str) -> Result<String, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("spacer sequence cannot be empty".to_string());
    }
    if trimmed.contains('U') || trimmed.contains('u') {
        return Err("spacer must use DNA alphabet (use 'T' instead of 'U')".to_string());
    }
    let upper = trimmed.to_ascii_uppercase();
    if let Some(ch) = upper.chars().find(|c| !matches!(c, 'A' | 'C' | 'G' | 'T' | 'N')) {
        return Err(format!("contains invalid character '{ch}'; expected A, C, G, T"));
    }
    Ok(upper)
}

/// Validate PAM motif: non-empty, IUPAC DNA characters only.
pub fn validate_manifest_pam(raw: &str) -> Result<String, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("PAM motif cannot be empty".to_string());
    }
    if trimmed.contains('U') || trimmed.contains('u') {
        return Err("PAM must use DNA alphabet (use 'T' instead of 'U')".to_string());
    }
    let upper = trimmed.to_ascii_uppercase();
    if let Some(ch) = upper.chars().find(|c| {
        !matches!(
            c,
            'A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N'
        )
    }) {
        return Err(format!("contains invalid character '{ch}'; expected IUPAC DNA code"));
    }
    Ok(upper)
}

impl BatchManifest {
    /// Parse batch manifest from raw TSV content string.
    ///
    /// If `base_dir` is provided, relative paths for `r1`, `r2`, and `intended_tsv`
    /// will be resolved relative to `base_dir`.
    pub fn parse_str(content: &str, base_dir: Option<&Path>) -> Result<Self, BatchManifestError> {
        let lines: Vec<&str> = content
            .lines()
            .map(|l| l.trim_end_matches('\r'))
            .filter(|l| !l.trim().is_empty())
            .collect();

        if lines.is_empty() {
            return Err(BatchManifestError::EmptyManifest);
        }

        let header_line = lines[0];
        let headers: Vec<&str> = header_line.split('\t').map(|s| s.trim()).collect();
        let expected_cols = headers.len();

        let clone_id_idx = headers.iter().position(|&h| h == "clone_id");
        let r1_idx = headers.iter().position(|&h| h == "r1");
        let r2_idx = headers.iter().position(|&h| h == "r2");
        let intended_idx = headers.iter().position(|&h| h == "intended_tsv");
        let spacer_idx = headers.iter().position(|&h| h == "spacer");
        let pam_idx = headers.iter().position(|&h| h == "pam");

        let mut missing = Vec::new();
        if clone_id_idx.is_none() {
            missing.push("clone_id");
        }
        if r1_idx.is_none() {
            missing.push("r1");
        }
        if !missing.is_empty() {
            return Err(BatchManifestError::MissingRequiredColumns(missing.join(", ")));
        }

        let clone_id_idx = clone_id_idx.unwrap();
        let r1_idx = r1_idx.unwrap();

        let mut entries = Vec::new();
        let mut seen_ids = HashSet::new();

        for (idx, line) in lines.iter().enumerate().skip(1) {
            let line_num = idx + 1;
            let parts: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();
            if parts.is_empty() || (parts.len() == 1 && parts[0].is_empty()) {
                continue;
            }

            // Strictly enforce row column width
            if parts.len() != expected_cols {
                return Err(BatchManifestError::RaggedRow {
                    line: line_num,
                    expected: expected_cols,
                    actual: parts.len(),
                });
            }

            let clone_id_raw = parts.get(clone_id_idx).copied().unwrap_or("");
            if clone_id_raw.is_empty() {
                return Err(BatchManifestError::EmptyCloneId { line: line_num });
            }

            let clone_id = validate_clone_id(clone_id_raw).map_err(|reason| {
                BatchManifestError::InvalidCloneId {
                    line: line_num,
                    clone_id: clone_id_raw.to_string(),
                    reason,
                }
            })?;

            if !seen_ids.insert(clone_id.clone()) {
                return Err(BatchManifestError::DuplicateCloneId {
                    line: line_num,
                    clone_id,
                });
            }

            let r1_str = parts.get(r1_idx).copied().unwrap_or("");
            if r1_str.is_empty() {
                return Err(BatchManifestError::EmptyR1 {
                    line: line_num,
                    clone_id,
                });
            }

            let resolve_path = |raw: &str| -> PathBuf {
                let p = PathBuf::from(raw);
                if p.is_relative() {
                    if let Some(base) = base_dir {
                        return base.join(p);
                    }
                }
                p
            };

            let r1 = resolve_path(r1_str);

            let r2 = r2_idx
                .and_then(|i| parts.get(i).copied())
                .and_then(|s| {
                    if s.is_empty() || s == "." || s.eq_ignore_ascii_case("na") {
                        None
                    } else {
                        Some(resolve_path(s))
                    }
                });

            let intended_tsv = intended_idx
                .and_then(|i| parts.get(i).copied())
                .and_then(|s| {
                    if s.is_empty() || s == "." || s.eq_ignore_ascii_case("na") {
                        None
                    } else {
                        Some(resolve_path(s))
                    }
                });

            let spacer = match spacer_idx.and_then(|i| parts.get(i).copied()) {
                Some(s) if !s.is_empty() && s != "." && !s.eq_ignore_ascii_case("na") => {
                    let validated = validate_manifest_spacer(s).map_err(|reason| {
                        BatchManifestError::InvalidSpacer {
                            line: line_num,
                            clone_id: clone_id.clone(),
                            reason,
                        }
                    })?;
                    Some(validated)
                }
                _ => None,
            };

            let pam = match pam_idx.and_then(|i| parts.get(i).copied()) {
                Some(s) if !s.is_empty() && s != "." && !s.eq_ignore_ascii_case("na") => {
                    let validated = validate_manifest_pam(s).map_err(|reason| {
                        BatchManifestError::InvalidPam {
                            line: line_num,
                            clone_id: clone_id.clone(),
                            reason,
                        }
                    })?;
                    Some(validated)
                }
                _ => None,
            };

            entries.push(BatchManifestEntry {
                clone_id,
                r1,
                r2,
                intended_tsv,
                spacer,
                pam,
            });
        }

        if entries.is_empty() {
            return Err(BatchManifestError::EmptyManifest);
        }

        Ok(Self { entries })
    }

    /// Preflight check that all referenced input files exist on disk.
    pub fn preflight_files(&self) -> Result<(), BatchManifestError> {
        for entry in &self.entries {
            if !entry.r1.exists() {
                return Err(BatchManifestError::InputFileNotFound {
                    clone_id: entry.clone_id.clone(),
                    path: entry.r1.clone(),
                    field: "r1",
                });
            }
            if let Some(ref r2) = entry.r2 {
                if !r2.exists() {
                    return Err(BatchManifestError::InputFileNotFound {
                        clone_id: entry.clone_id.clone(),
                        path: r2.clone(),
                        field: "r2",
                    });
                }
            }
            if let Some(ref it) = entry.intended_tsv {
                if !it.exists() {
                    return Err(BatchManifestError::InputFileNotFound {
                        clone_id: entry.clone_id.clone(),
                        path: it.clone(),
                        field: "intended_tsv",
                    });
                }
            }
        }
        Ok(())
    }

    /// Read, parse, and preflight batch manifest from a file path.
    pub fn parse_file(path: impl AsRef<Path>) -> Result<Self, BatchManifestError> {
        let path = path.as_ref();
        if !path.exists() {
            return Err(BatchManifestError::FileNotFound(path.to_path_buf()));
        }
        let content = fs::read_to_string(path)
            .map_err(|e| BatchManifestError::Io(format!("{}: {e}", path.display())))?;
        let base_dir = path.parent();
        let manifest = Self::parse_str(&content, base_dir)?;
        manifest.preflight_files()?;
        Ok(manifest)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_valid_manifest_mixed_pe_se() {
        let content = "clone_id\tr1\tr2\tintended_tsv\tspacer\tpam\n\
                       C01\t/data/c1_1.fq\t/data/c1_2.fq\t/data/intended.tsv\tATCGATCGATCGATCGATCG\tNGG\n\
                       C02\t/data/c2_1.fq\t.\t.\t.\t.\n\
                       C03\t/data/c3_1.fq\tNA\tNA\tATCGATCGATCGATCGATCG\tTTTV\n";

        let manifest = BatchManifest::parse_str(content, None).unwrap();
        assert_eq!(manifest.entries.len(), 3);

        let e1 = &manifest.entries[0];
        assert_eq!(e1.clone_id, "C01");
        assert_eq!(e1.r1, PathBuf::from("/data/c1_1.fq"));
        assert_eq!(e1.r2, Some(PathBuf::from("/data/c1_2.fq")));
        assert_eq!(e1.intended_tsv, Some(PathBuf::from("/data/intended.tsv")));
        assert_eq!(e1.spacer, Some("ATCGATCGATCGATCGATCG".to_string()));
        assert_eq!(e1.pam, Some("NGG".to_string()));

        let e2 = &manifest.entries[1];
        assert_eq!(e2.clone_id, "C02");
        assert_eq!(e2.r1, PathBuf::from("/data/c2_1.fq"));
        assert_eq!(e2.r2, None);
        assert_eq!(e2.intended_tsv, None);
        assert_eq!(e2.spacer, None);
        assert_eq!(e2.pam, None);
    }

    #[test]
    fn test_reject_ragged_row() {
        let content = "clone_id\tr1\tr2\n\
                       C01\tr1.fq\n"; // 2 columns instead of 3
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert!(matches!(
            err,
            BatchManifestError::RaggedRow {
                line: 2,
                expected: 3,
                actual: 2
            }
        ));
    }

    #[test]
    fn test_reject_invalid_clone_ids() {
        // Path traversal
        let content = "clone_id\tr1\n../escape\tr1.fq\n";
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert!(matches!(err, BatchManifestError::InvalidCloneId { .. }));

        // Slashing
        let content = "clone_id\tr1\nsub/dir\tr1.fq\n";
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert!(matches!(err, BatchManifestError::InvalidCloneId { .. }));

        // Dot
        let content = "clone_id\tr1\n.\tr1.fq\n";
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert!(matches!(err, BatchManifestError::InvalidCloneId { .. }));
    }

    #[test]
    fn test_reject_invalid_spacer_and_pam() {
        // Spacer with U
        let content = "clone_id\tr1\tspacer\nC01\tr1.fq\tAUCGATCG\n";
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert!(matches!(err, BatchManifestError::InvalidSpacer { .. }));

        // PAM with invalid character
        let content = "clone_id\tr1\tpam\nC01\tr1.fq\t@@@\n";
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert!(matches!(err, BatchManifestError::InvalidPam { .. }));
    }

    #[test]
    fn test_reject_duplicate_clone_id() {
        let content = "clone_id\tr1\nC01\tr1.fq\nC01\tr2.fq\n";
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert!(matches!(
            err,
            BatchManifestError::DuplicateCloneId { clone_id, .. } if clone_id == "C01"
        ));
    }

    #[test]
    fn test_reject_empty_manifest() {
        let content = "clone_id\tr1\n";
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert_eq!(err, BatchManifestError::EmptyManifest);
    }

    #[test]
    fn test_reject_missing_required_headers() {
        let content = "sample\treads\nS1\tr.fq\n";
        let err = BatchManifest::parse_str(content, None).unwrap_err();
        assert!(matches!(
            err,
            BatchManifestError::MissingRequiredColumns(_)
        ));
    }
}
