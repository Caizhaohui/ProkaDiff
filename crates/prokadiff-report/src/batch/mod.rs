//! Multi-clone batch audit mode domain models, parsing, recurrence analysis, and reporting.

pub mod manifest;
pub mod recurrence;
pub mod render;

pub use manifest::{BatchManifest, BatchManifestEntry, BatchManifestError};
pub use recurrence::{
    align_batch_variants, BatchVariantRecord, CloneVariantObservation, RecurrenceStatus,
    VariantLocusKey,
};
pub use render::{
    extract_clone_observations, render_batch_report_md, render_batch_variants_tsv,
    write_batch_report_md, write_batch_variants_tsv, BatchAuditSummary, CloneBatchSummary,
};
