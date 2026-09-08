use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::error::{EvidenceError, Result};

use super::FastaRecord;

pub(super) fn validate_reference_ids(
    records: &[FastaRecord],
    source: &Path,
    ids: &mut HashMap<String, PathBuf>,
) -> Result<()> {
    for record in records {
        if let Some(first_path) = ids.get(&record.name) {
            return Err(EvidenceError::DuplicateReferenceId {
                id: record.name.clone(),
                first_path: first_path.clone(),
                duplicate_path: source.to_path_buf(),
            });
        }
        ids.insert(record.name.clone(), source.to_path_buf());
    }
    Ok(())
}
