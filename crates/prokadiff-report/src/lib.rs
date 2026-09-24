mod markdown;
mod schema;
mod tables;

pub use markdown::write_markdown_report;
pub use schema::SchemaVersion;
pub use tables::{write_edit_outcomes_tsv, write_post_edit_variants_tsv, write_provenance_tsv};

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use prokadiff_classify::{
    write_offtarget_sites_tsv, AuditResult, ClassifiedMutation, ClassifyResult, EventDisplay,
    MutationClass, RefContig,
};
use prokadiff_gd::{GdEntry, GdKind};

/// Write classified post-subtract mutations. `hypothesis` column omitted when `include_hypothesis` is false.
///
/// Column order: existing product columns through `distance_to_site`, then `side2_seq_id` /
/// `side2_position` (empty except JC), then optional `hypothesis`.
pub fn write_unintended_tsv(
    path: impl AsRef<Path>,
    rows: &[ClassifiedMutation],
    editor: &str,
    include_hypothesis: bool,
    refs: &[RefContig],
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write!(
        w,
        "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position"
    )?;
    if include_hypothesis {
        write!(w, "\thypothesis")?;
    }
    writeln!(w)?;
    for row in rows {
        let (seq_id, pos, end) = coords(&row.entry);
        let (ref_a, alt) = alleles(&row.entry, refs);
        let (side2_id, side2_pos) = side2(&row.entry);
        write!(
            w,
            "{seq_id}\t{pos}\t{end}\t{}\t{ref_a}\t{alt}\t{}\t{editor}\t{}\t{}\t{}\t{side2_id}\t{side2_pos}",
            row.entry.kind.as_str(),
            row.class.as_str(),
            row.pam_profile.as_deref().unwrap_or(""),
            opt_u32(row.offtarget_mismatch),
            opt_u64(row.distance_to_site),
        )?;
        if include_hypothesis {
            write!(w, "\t{}", row.hypothesis.as_deref().unwrap_or(""))?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}

fn write_unintended_from_audit(
    path: impl AsRef<Path>,
    audit: &AuditResult,
    schema_version: SchemaVersion,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write!(
        w,
        "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position"
    )?;
    if audit.hypothesis_enabled {
        write!(w, "\thypothesis")?;
    }
    if schema_version.is_v2() {
        write!(w, "\tevent_id")?;
    }
    writeln!(w)?;
    for row in &audit.unintended {
        let (seq_id, pos, end) = coords(&row.entry);
        let display = event_display(audit, row.event_id.as_ref())?;
        let (side2_id, side2_pos) = side2(&row.entry);
        write!(
            w,
            "{seq_id}\t{pos}\t{end}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{side2_id}\t{side2_pos}",
            row.entry.kind.as_str(),
            display.unintended_reference,
            display.unintended_alternate,
            row.class.as_str(),
            audit.sample.editor,
            row.pam_profile.as_deref().unwrap_or(""),
            opt_u32(row.offtarget_mismatch),
            opt_u64(row.distance_to_site),
        )?;
        if audit.hypothesis_enabled {
            write!(w, "\t{}", row.hypothesis.as_deref().unwrap_or(""))?;
        }
        if schema_version.is_v2() {
            write!(
                w,
                "\t{}",
                row.event_id
                    .as_ref()
                    .map_or("NA", |event_id| event_id.as_str())
            )?;
        }
        writeln!(w)?;
    }
    w.flush()
}

/// Write the run summary.  Includes:
///  - FIX-015 edit-level fields: `intended_edits_complete`, `intended_edits_partial`,
///    `intended_edits_missing`, `intended_events_observed`
///  - Deprecated (backward-compat): `intended_declared`, `intended_observed`,
///    `intended_status`, `intended_missing`
///  - FIX-018 provenance: validation status for scoring components
pub fn write_summary(
    path: impl AsRef<Path>,
    result: &ClassifyResult,
    intended_provided: bool,
    editor: &str,
) -> std::io::Result<()> {
    let mut n_s = 0usize;
    let mut n_n = 0usize;
    let mut n_c = 0usize;
    for r in &result.unintended {
        match r.class {
            MutationClass::Structural => n_s += 1,
            MutationClass::NearHomolog => n_n += 1,
            MutationClass::ScatteredSnv => n_c += 1,
        }
    }
    let (declared, observed, status, missing) = intended_summary_fields(intended_provided, result);
    // FIX-015: edit-level breakdown (requires assess_intended_edits call upstream)
    let (edits_complete, edits_partial, edits_missing, events_observed) =
        intended_edit_level_fields(intended_provided, result);
    let text = format!(
        "editor\t{editor}\n\
# -- Intended edit summary (FIX-015) --\n\
intended_provided\t{}\n\
intended_edits_declared\t{declared}\n\
intended_edits_complete\t{edits_complete}\n\
intended_edits_partial\t{edits_partial}\n\
intended_edits_missing\t{edits_missing}\n\
intended_events_observed\t{events_observed}\n\
# -- Deprecated fields (backward compatibility) --\n\
intended_declared\t{declared}\n\
intended_observed\t{observed}\n\
intended_status\t{status}\n\
intended_missing\t{missing}\n\
# -- Mutation class counts --\n\
structural\t{n_s}\n\
near_homolog\t{n_n}\n\
scattered_snv\t{n_c}\n\
starter_vs_ref_mutations\t{}\n\
# -- FIX-018 validation status --\n\
offtarget_search_validation_status\tvalidated_via_self_test\n\
cfd_validation_status\tdisabled\n\
hsu_validation_status\texperimental\n\
bulge_validation_status\texperimental\n",
        if intended_provided { "yes" } else { "no" },
        result.starter_vs_ref,
    );
    std::fs::write(path, text)?;
    Ok(())
}

fn write_summary_from_audit(path: impl AsRef<Path>, audit: &AuditResult) -> std::io::Result<()> {
    let summary = &audit.summary;
    let display =
        |value: Option<usize>| value.map_or_else(|| "NA".to_string(), |value| value.to_string());
    let declared = display(summary.intended_declared);
    let observed = display(summary.intended_events_observed);
    let status = summary
        .intended_status
        .map_or("NA", |status| status.as_str());
    let missing = display(summary.intended_missing);
    let edits_complete = display(summary.intended_edits_complete);
    let edits_partial = display(summary.intended_edits_partial);
    let edits_missing = display(summary.intended_edits_missing);
    let events = display(summary.intended_events_observed);
    let text = format!(
        "editor\t{}\n\
# -- Intended edit summary (FIX-015) --\n\
intended_provided\t{}\n\
intended_edits_declared\t{declared}\n\
intended_edits_complete\t{edits_complete}\n\
intended_edits_partial\t{edits_partial}\n\
intended_edits_missing\t{edits_missing}\n\
intended_events_observed\t{events}\n\
# -- Deprecated fields (backward compatibility) --\n\
intended_declared\t{declared}\n\
intended_observed\t{observed}\n\
intended_status\t{status}\n\
intended_missing\t{missing}\n\
# -- Mutation class counts --\n\
structural\t{}\n\
near_homolog\t{}\n\
scattered_snv\t{}\n\
starter_vs_ref_mutations\t{}\n\
# -- FIX-018 validation status --\n\
offtarget_search_validation_status\tvalidated_via_self_test\n\
cfd_validation_status\tdisabled\n\
hsu_validation_status\texperimental\n\
bulge_validation_status\texperimental\n",
        audit.sample.editor,
        if summary.intended_provided {
            "yes"
        } else {
            "no"
        },
        summary.structural_count,
        summary.near_homolog_count,
        summary.scattered_snv_count,
        summary.starter_vs_reference,
    );
    std::fs::write(path, text)
}

const PRODUCT_OUTPUT_NAMES: [&str; 8] = [
    "unintended.tsv",
    "summary.txt",
    "offtarget_sites.tsv",
    "mutation_offtarget_links.tsv",
    "report.md",
    "edit_outcomes.tsv",
    "post_edit_variants.tsv",
    "provenance.tsv",
];

#[derive(Debug, thiserror::Error)]
pub enum ProductOutputError {
    #[error("product output I/O failed at {}: {source}", path.display())]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error(
        "product publication failed at {}; rollback failed at {}; backup retained at {}; staging retained at {}: {source}",
        publish_path.display(),
        recovery_path.display(),
        backup_path.display(),
        staging_path.display(),
    )]
    Rollback {
        publish_path: PathBuf,
        recovery_path: PathBuf,
        backup_path: PathBuf,
        staging_path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error(
        "product publication cleanup failed while {operation} at {}; backup retained at {}; staging retained at {}: {source}",
        path.display(),
        backup_path.display(),
        staging_path.display(),
    )]
    Cleanup {
        operation: &'static str,
        path: PathBuf,
        backup_path: PathBuf,
        staging_path: PathBuf,
        #[source]
        source: std::io::Error,
    },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PublicationOperation {
    Publish,
    Restore,
}

#[cfg(test)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PublicationFault {
    PublishAt(usize),
    PublishAtThenRestoreAt { publish: usize, restore: usize },
}

#[cfg(test)]
std::thread_local! {
    static PUBLICATION_FAULT: std::cell::Cell<Option<PublicationFault>> = const {
        std::cell::Cell::new(None)
    };
}

#[cfg(test)]
fn set_publication_fault(fault: Option<PublicationFault>) {
    PUBLICATION_FAULT.with(|configured| configured.set(fault));
}

fn inject_publication_failure(
    _operation: PublicationOperation,
    _index: usize,
    _path: &Path,
) -> std::io::Result<()> {
    #[cfg(test)]
    {
        let matches = PUBLICATION_FAULT.with(|configured| match configured.get() {
            Some(PublicationFault::PublishAt(index)) => {
                _operation == PublicationOperation::Publish && _index == index
            }
            Some(PublicationFault::PublishAtThenRestoreAt { publish, restore }) => {
                (_operation == PublicationOperation::Publish && _index == publish)
                    || (_operation == PublicationOperation::Restore && _index == restore)
            }
            None => false,
        });
        if matches {
            return Err(std::io::Error::other(format!(
                "injected {:?} failure at {}",
                _operation,
                _path.display()
            )));
        }
    }
    Ok(())
}

fn output_io(path: impl Into<PathBuf>, source: std::io::Error) -> ProductOutputError {
    ProductOutputError::Io {
        path: path.into(),
        source,
    }
}

fn rollback_error(
    publish_path: &Path,
    recovery_path: impl Into<PathBuf>,
    backup_path: &Path,
    staging_path: &Path,
    source: std::io::Error,
) -> ProductOutputError {
    ProductOutputError::Rollback {
        publish_path: publish_path.to_path_buf(),
        recovery_path: recovery_path.into(),
        backup_path: backup_path.to_path_buf(),
        staging_path: staging_path.to_path_buf(),
        source,
    }
}

fn cleanup_error(
    operation: &'static str,
    path: impl Into<PathBuf>,
    backup_path: &Path,
    staging_path: &Path,
    source: std::io::Error,
) -> ProductOutputError {
    ProductOutputError::Cleanup {
        operation,
        path: path.into(),
        backup_path: backup_path.to_path_buf(),
        staging_path: staging_path.to_path_buf(),
        source,
    }
}

fn restore_previous_product_set(
    outdir: &Path,
    staging: &Path,
    backup: &Path,
    backed_up: &[&str],
    published: &[&str],
    publish_path: &Path,
) -> Result<(), ProductOutputError> {
    for name in published.iter().rev() {
        let path = outdir.join(name);
        if let Err(source) = fs::remove_file(&path) {
            return Err(rollback_error(publish_path, path, backup, staging, source));
        }
    }

    for (index, name) in backed_up.iter().enumerate() {
        let source_path = backup.join(name);
        let destination = outdir.join(name);
        if let Err(source) =
            inject_publication_failure(PublicationOperation::Restore, index, &destination)
        {
            return Err(rollback_error(
                publish_path,
                destination,
                backup,
                staging,
                source,
            ));
        }
        if let Err(source) = fs::copy(&source_path, &destination) {
            return Err(rollback_error(
                publish_path,
                destination,
                backup,
                staging,
                source,
            ));
        }
    }

    for name in backed_up {
        let source_path = backup.join(name);
        let destination = outdir.join(name);
        let source_size = match fs::metadata(&source_path) {
            Ok(metadata) => metadata.len(),
            Err(source) => {
                return Err(rollback_error(
                    publish_path,
                    source_path,
                    backup,
                    staging,
                    source,
                ));
            }
        };
        let destination_size = match fs::metadata(&destination) {
            Ok(metadata) if metadata.is_file() => metadata.len(),
            Ok(_) => {
                return Err(rollback_error(
                    publish_path,
                    destination,
                    backup,
                    staging,
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "restored product output is not a regular file",
                    ),
                ));
            }
            Err(source) => {
                return Err(rollback_error(
                    publish_path,
                    destination,
                    backup,
                    staging,
                    source,
                ));
            }
        };
        if source_size != destination_size {
            return Err(rollback_error(
                publish_path,
                destination,
                backup,
                staging,
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "restored product output size does not match backup",
                ),
            ));
        }
    }
    Ok(())
}

pub fn write_product_outputs(
    audit: &AuditResult,
    schema_version: SchemaVersion,
    outdir: impl AsRef<Path>,
) -> Result<(), ProductOutputError> {
    let outdir = outdir.as_ref();
    fs::create_dir_all(outdir).map_err(|source| output_io(outdir, source))?;
    let staging = outdir.join(format!(
        ".prokadiff-product-stage-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos()
    ));
    fs::create_dir(&staging).map_err(|source| output_io(&staging, source))?;
    write_product_set(&staging, audit, schema_version)
        .map_err(|source| output_io(&staging, source))?;
    if let Some(path) = PRODUCT_OUTPUT_NAMES
        .iter()
        .map(|name| outdir.join(name))
        .find(|path| path.is_dir())
    {
        return Err(output_io(
            path.clone(),
            std::io::Error::new(
                std::io::ErrorKind::AlreadyExists,
                format!("product output path is a directory: {}", path.display()),
            ),
        ));
    }
    let backup = outdir.join(format!(
        ".prokadiff-product-backup-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos()
    ));
    fs::create_dir(&backup).map_err(|source| output_io(&backup, source))?;
    let mut backed_up = Vec::new();
    for name in PRODUCT_OUTPUT_NAMES {
        let destination = outdir.join(name);
        if destination.exists() {
            if let Err(error) = fs::rename(&destination, backup.join(name)) {
                restore_previous_product_set(
                    outdir,
                    &staging,
                    &backup,
                    &backed_up,
                    &[],
                    &destination,
                )?;
                if let Err(source) = fs::remove_dir_all(&staging) {
                    return Err(cleanup_error(
                        "removing staging after backup rollback",
                        &staging,
                        &backup,
                        &staging,
                        source,
                    ));
                }
                if let Err(source) = fs::remove_dir_all(&backup) {
                    return Err(cleanup_error(
                        "removing backup after backup rollback",
                        &backup,
                        &backup,
                        &staging,
                        source,
                    ));
                }
                return Err(output_io(destination, error));
            }
            backed_up.push(name);
        }
    }
    let mut published = Vec::new();
    for (index, name) in PRODUCT_OUTPUT_NAMES.iter().enumerate() {
        let destination = outdir.join(name);
        let publish =
            inject_publication_failure(PublicationOperation::Publish, index, &destination)
                .and_then(|()| fs::rename(staging.join(name), &destination));
        if let Err(error) = publish {
            restore_previous_product_set(
                outdir,
                &staging,
                &backup,
                &backed_up,
                &published,
                &destination,
            )?;
            if let Err(source) = fs::remove_dir_all(&staging) {
                return Err(cleanup_error(
                    "removing staging after publication rollback",
                    &staging,
                    &backup,
                    &staging,
                    source,
                ));
            }
            if let Err(source) = fs::remove_dir_all(&backup) {
                return Err(cleanup_error(
                    "removing backup after verified publication rollback",
                    &backup,
                    &backup,
                    &staging,
                    source,
                ));
            }
            return Err(output_io(destination, error));
        }
        published.push(name);
    }
    if let Err(source) = fs::remove_dir(&staging) {
        return Err(cleanup_error(
            "removing empty staging after publication",
            &staging,
            &backup,
            &staging,
            source,
        ));
    }
    if let Err(source) = fs::remove_dir_all(&backup) {
        return Err(cleanup_error(
            "removing backup after publication",
            &backup,
            &backup,
            &staging,
            source,
        ));
    }
    Ok(())
}

fn write_product_set(
    outdir: &Path,
    audit: &AuditResult,
    schema_version: SchemaVersion,
) -> std::io::Result<()> {
    write_unintended_from_audit(outdir.join("unintended.tsv"), audit, schema_version)?;
    write_summary_from_audit(outdir.join("summary.txt"), audit)?;
    write_offtarget_sites_tsv(&audit.guide_sites, outdir.join("offtarget_sites.tsv"))?;
    tables::write_mutation_offtarget_links_from_audit(
        outdir.join("mutation_offtarget_links.tsv"),
        audit,
        schema_version,
    )?;
    write_markdown_report(outdir.join("report.md"), audit, schema_version)?;
    tables::write_edit_outcomes_from_audit(
        outdir.join("edit_outcomes.tsv"),
        audit,
        schema_version,
    )?;
    tables::write_post_edit_variants_from_audit(
        outdir.join("post_edit_variants.tsv"),
        audit,
        schema_version,
    )?;
    write_provenance_tsv(outdir.join("provenance.tsv"), &audit.provenance)
}

pub(crate) fn event_display<'a>(
    audit: &'a AuditResult,
    event_id: Option<&prokadiff_classify::EventId>,
) -> std::io::Result<&'a EventDisplay> {
    let event_id = event_id.ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "missing EventId for product row",
        )
    })?;
    audit.event_display.get(event_id).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("missing event display for {event_id}"),
        )
    })
}

/// Returns (complete, partial, missing, events_observed) for edit-level summary.
///
/// Uses `ClassifyResult.intended_edit_assessments` when available;
/// falls back to approximate event-count logic for backward compatibility.
fn intended_edit_level_fields(
    provided: bool,
    result: &ClassifyResult,
) -> (String, String, String, String) {
    if !provided {
        return ("NA".into(), "NA".into(), "NA".into(), "NA".into());
    }
    if let Some(assessments) = &result.intended_edit_assessments {
        use prokadiff_classify::IntendedEditStatus;
        let complete = assessments
            .iter()
            .filter(|a| a.status == IntendedEditStatus::Complete)
            .count();
        let partial = assessments
            .iter()
            .filter(|a| a.status == IntendedEditStatus::Partial)
            .count();
        let missing = assessments
            .iter()
            .filter(|a| a.status == IntendedEditStatus::Missing)
            .count();
        let events = result.intended_event_ids.len();
        (
            complete.to_string(),
            partial.to_string(),
            missing.to_string(),
            events.to_string(),
        )
    } else {
        // Approximate from event count (pre-FIX-015 path)
        let declared = result.intended_declared;
        let observed = result.intended_observed.len();
        if declared == 0 {
            return ("NA".into(), "NA".into(), "NA".into(), "NA".into());
        }
        // Treat observed events as complete edits (approximate)
        (
            observed.to_string(),
            "0".into(),
            declared.saturating_sub(observed).to_string(),
            observed.to_string(),
        )
    }
}

fn intended_summary_fields(
    provided: bool,
    result: &ClassifyResult,
) -> (String, String, &'static str, String) {
    if !provided {
        return ("NA".into(), "NA".into(), "NA", "NA".into());
    }
    let declared = result.intended_declared;
    let observed = result.intended_observed.len();
    if declared == 0 {
        return (
            declared.to_string(),
            observed.to_string(),
            "NA",
            "NA".into(),
        );
    }
    // FIX-015: status is now determined by edit-level assessment, not raw event counts
    // Use assessments if available; fall back to event count for backward compat.
    let status = if let Some(assessments) = &result.intended_edit_assessments {
        use prokadiff_classify::IntendedEditStatus;
        let complete = assessments
            .iter()
            .filter(|a| a.status == IntendedEditStatus::Complete)
            .count();
        let missing = assessments
            .iter()
            .filter(|a| a.status == IntendedEditStatus::Missing)
            .count();
        if complete == declared {
            "all_observed"
        } else if missing == declared {
            "none_observed"
        } else {
            "partial"
        }
    } else {
        // Legacy event-count comparison
        if observed >= declared {
            "all_observed"
        } else if observed == 0 {
            "none_observed"
        } else {
            "partial"
        }
    };
    (
        declared.to_string(),
        observed.to_string(),
        status,
        declared.saturating_sub(observed).to_string(),
    )
}

fn opt_u32(v: Option<u32>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}

fn opt_u64(v: Option<u64>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}

fn coords(e: &GdEntry) -> (String, u64, u64) {
    let seq = e.seq_id().unwrap_or("").to_string();
    let pos = e.position().unwrap_or(0);
    if e.kind == GdKind::Del || e.kind == GdKind::Sub || e.kind == GdKind::Inv {
        let size = e
            .fields
            .get(2)
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(1)
            .max(1);
        (seq, pos, pos.saturating_add(size - 1))
    } else {
        (seq, pos, pos)
    }
}

fn alleles<'a>(e: &'a GdEntry, refs: &[RefContig]) -> (String, &'a str) {
    match e.kind {
        GdKind::Snp => {
            let seq_id = e.seq_id().unwrap_or("");
            let pos = e.position().unwrap_or(0);
            (
                snp_ref_base(seq_id, pos, refs),
                e.fields.get(2).map(String::as_str).unwrap_or("."),
            )
        }
        GdKind::Sub => {
            let seq_id = e.seq_id().unwrap_or("");
            let pos = e.position().unwrap_or(0);
            let size = e
                .fields
                .get(2)
                .and_then(|s| s.parse::<usize>().ok())
                .unwrap_or(1);
            (
                sub_ref_seq(seq_id, pos, size, refs),
                e.fields.get(3).map(String::as_str).unwrap_or("."),
            )
        }
        GdKind::Ins => (
            ".".into(),
            e.fields.get(2).map(String::as_str).unwrap_or("."),
        ),
        _ => (".".into(), "."),
    }
}

fn sub_ref_seq(seq_id: &str, pos: u64, size: usize, refs: &[RefContig]) -> String {
    if pos == 0 {
        return ".".into();
    }
    let start_idx = (pos - 1) as usize;
    for c in refs {
        if c.name == seq_id && start_idx < c.seq.len() {
            let end_idx = (start_idx + size).min(c.seq.len());
            let slice = &c.seq[start_idx..end_idx];
            return slice
                .iter()
                .map(|&b| (b as char).to_ascii_uppercase())
                .collect();
        }
    }
    ".".into()
}

fn snp_ref_base(seq_id: &str, pos: u64, refs: &[RefContig]) -> String {
    if pos == 0 {
        return ".".into();
    }
    let idx = (pos - 1) as usize;
    for c in refs {
        if c.name == seq_id {
            if let Some(&b) = c.seq.get(idx) {
                return (b as char).to_ascii_uppercase().to_string();
            }
        }
    }
    ".".into()
}

fn side2(e: &GdEntry) -> (String, String) {
    if e.kind == GdKind::Jc {
        (
            e.fields.get(3).cloned().unwrap_or_default(),
            e.fields.get(4).cloned().unwrap_or_default(),
        )
    } else {
        (String::new(), String::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use prokadiff_classify::{
        build_complete_audit_result, AnalysisProvenance, CandidateSearchStatus, ClassifiedMutation,
        MutationClass, RefContig, SampleMetadata,
    };
    use prokadiff_gd::GdEntry;
    use std::collections::BTreeMap;

    const TSV_HEAD_NO_HYP: &str = "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position";
    const TSV_HEAD_WITH_HYP: &str = "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position\thypothesis";

    fn scratch(name: &str) -> std::path::PathBuf {
        let dir =
            std::env::temp_dir().join(format!("prokdiff-report-{}-{}", name, std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn snp_mut(pos: u64, alt: &str, class: MutationClass) -> ClassifiedMutation {
        ClassifiedMutation {
            entry: GdEntry::snp(1, "chr", pos, alt),
            class,
            pam_profile: None,
            offtarget_mismatch: None,
            distance_to_site: None,
            hypothesis: Some("sos_widney2014".into()),
            event_id: None,
        }
    }

    fn classify_result(
        declared: usize,
        observed: usize,
        unintended: Vec<ClassifiedMutation>,
        starter_vs_ref: usize,
    ) -> ClassifyResult {
        ClassifyResult {
            unintended,
            intended_observed: (0..observed)
                .map(|i| GdEntry::snp(i as u32 + 10, "chr", 1000 + i as u64, "A"))
                .collect(),
            intended_event_ids: Vec::new(),
            intended_declared: declared,
            starter_vs_ref,
            intended_edit_assessments: None,
            differential_events: Vec::new(),
            candidate_search: Default::default(),
            associations: Vec::new(),
        }
    }

    fn kv(text: &str, key: &str) -> String {
        for line in text.lines() {
            let mut parts = line.splitn(2, '\t');
            if parts.next() == Some(key) {
                return parts.next().unwrap_or("").to_string();
            }
        }
        panic!("missing summary key {key} in:\n{text}");
    }

    fn data_cols(path: &std::path::Path) -> Vec<String> {
        let text = std::fs::read_to_string(path).unwrap();
        text.lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .map(str::to_string)
            .collect()
    }

    fn publication_audit(editor: &str) -> AuditResult {
        build_complete_audit_result(
            SampleMetadata {
                starter_names: vec!["starter.fq".into()],
                edited_names: vec!["edited.fq".into()],
                reference_names: vec!["reference.fa".into()],
                editor: editor.into(),
                spacer: None,
                pam: None,
                threads: 1,
            },
            false,
            Vec::new(),
            Vec::new(),
            0,
            false,
            &[],
            &[],
            &[],
            &[],
            CandidateSearchStatus::NotPerformed,
            AnalysisProvenance {
                prokadiff_version: "test".into(),
                git_commit: "test".into(),
                reference_sha256: None,
                bowtie2_version: None,
                offtarget_search_status: "NOT_REQUESTED".into(),
                cfd_scoring_status: "DISABLED".into(),
                hsu_scoring_status: "DISABLED".into(),
                bulge_search_status: "EXACT_UNGAPPED".into(),
                run_timestamp: "2026-09-24T00:00:00Z".into(),
            },
            &[],
            &[],
        )
        .expect("valid empty audit aggregate")
    }

    fn publication_scratch(name: &str) -> std::path::PathBuf {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("system time after epoch")
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "prokadiff-publication-{name}-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&dir).expect("create publication scratch directory");
        dir
    }

    fn product_bytes(outdir: &Path) -> BTreeMap<&'static str, Vec<u8>> {
        PRODUCT_OUTPUT_NAMES
            .iter()
            .map(|name| {
                (
                    *name,
                    std::fs::read(outdir.join(name)).expect("read complete product output"),
                )
            })
            .collect()
    }

    fn backup_paths(outdir: &Path) -> Vec<std::path::PathBuf> {
        std::fs::read_dir(outdir)
            .expect("read publication directory")
            .filter_map(Result::ok)
            .map(|entry| entry.path())
            .filter(|path| {
                path.file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.starts_with(".prokadiff-product-backup-"))
            })
            .collect()
    }

    struct PublicationFaultGuard;

    impl Drop for PublicationFaultGuard {
        fn drop(&mut self) {
            set_publication_fault(None);
        }
    }

    fn inject_fault(fault: PublicationFault) -> PublicationFaultGuard {
        set_publication_fault(Some(fault));
        PublicationFaultGuard
    }

    #[test]
    fn tsv_omits_hypothesis_column_when_disabled() {
        let path = scratch("no-hyp").join("unintended.tsv");
        let row = snp_mut(40, "C", MutationClass::ScatteredSnv);
        write_unintended_tsv(&path, &[row], "cas9", false, &[]).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(
            text.starts_with(&format!("{TSV_HEAD_NO_HYP}\n")),
            "header was: {:?}",
            text.lines().next()
        );
        assert!(!text.contains("hypothesis"));
        assert!(text.contains("scattered_snv"));
        assert!(text.contains("cas9"));
        assert!(text.contains("side2_seq_id"));
    }

    #[test]
    fn tsv_keeps_side2_columns_when_hypothesis_enabled() {
        let path = scratch("with-hyp").join("unintended.tsv");
        let row = snp_mut(40, "C", MutationClass::ScatteredSnv);
        write_unintended_tsv(&path, &[row], "cas9", true, &[]).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(text.lines().next().unwrap(), TSV_HEAD_WITH_HYP);
    }

    #[test]
    fn snp_tsv_ref_is_reference_base_at_1based_position() {
        let path = scratch("snp-ref").join("unintended.tsv");
        // 1-based pos 3 → 'G'
        let refs = [RefContig {
            name: "chr".into(),
            seq: b"ACGT".to_vec(),
        }];
        let row = snp_mut(3, "T", MutationClass::ScatteredSnv);
        write_unintended_tsv(&path, &[row], "cas9", false, &refs).unwrap();
        let data = std::fs::read_to_string(&path).unwrap();
        let body = data.lines().nth(1).unwrap();
        let cols: Vec<&str> = body.split('\t').collect();
        assert_eq!(cols[4], "G", "SNP ref column: {body}");
        assert_eq!(cols[5], "T");
        assert_eq!(cols[11], "", "SNP side2_seq_id must be empty");
        assert_eq!(cols[12], "", "SNP side2_position must be empty");
    }

    #[test]
    fn ins_tsv_ref_is_dot() {
        let path = scratch("ins-ref").join("unintended.tsv");
        let refs = [RefContig {
            name: "chr".into(),
            seq: b"ACGT".to_vec(),
        }];
        let row = ClassifiedMutation {
            entry: GdEntry::ins(1, "chr", 2, "AA"),
            class: MutationClass::ScatteredSnv,
            pam_profile: None,
            offtarget_mismatch: None,
            distance_to_site: None,
            hypothesis: None,
            event_id: None,
        };
        write_unintended_tsv(&path, &[row], "cas9", false, &refs).unwrap();
        let cols = data_cols(&path);
        assert_eq!(cols[4], ".");
        assert_eq!(cols[5], "AA");
    }

    #[test]
    fn jc_tsv_writes_side2_coordinates() {
        let path = scratch("jc-side2").join("unintended.tsv");
        let row = ClassifiedMutation {
            entry: GdEntry::jc(1, "chr", 100, "+", "plasmid", 55, "-", 0),
            class: MutationClass::Structural,
            pam_profile: None,
            offtarget_mismatch: None,
            distance_to_site: None,
            hypothesis: None,
            event_id: None,
        };
        write_unintended_tsv(&path, &[row], "cas9", false, &[]).unwrap();
        let cols = data_cols(&path);
        assert_eq!(cols[0], "chr");
        assert_eq!(cols[1], "100");
        assert_eq!(cols[3], "JC");
        assert_eq!(cols[11], "plasmid");
        assert_eq!(cols[12], "55");
    }

    #[test]
    fn summary_marks_intended_na_when_omitted() {
        let path = scratch("sum-na").join("summary.txt");
        let result = classify_result(
            0,
            0,
            vec![ClassifiedMutation {
                entry: GdEntry::snp(1, "chr", 40, "C"),
                class: MutationClass::NearHomolog,
                pam_profile: Some("NGG".into()),
                offtarget_mismatch: Some(0),
                distance_to_site: Some(7),
                hypothesis: None,
                event_id: None,
            }],
            50,
        );
        write_summary(&path, &result, false, "cas9").unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(kv(&text, "intended_provided"), "no");
        assert_eq!(kv(&text, "intended_declared"), "NA");
        assert_eq!(kv(&text, "intended_observed"), "NA");
        assert_eq!(kv(&text, "intended_status"), "NA");
        assert_eq!(kv(&text, "intended_missing"), "NA");
        assert_eq!(kv(&text, "near_homolog"), "1");
        assert_eq!(kv(&text, "starter_vs_ref_mutations"), "50");
    }

    #[test]
    fn summary_all_observed_when_declared_equals_observed() {
        let path = scratch("sum-all").join("summary.txt");
        let result = classify_result(2, 2, vec![], 0);
        write_summary(&path, &result, true, "cas9").unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(kv(&text, "intended_provided"), "yes");
        assert_eq!(kv(&text, "intended_declared"), "2");
        assert_eq!(kv(&text, "intended_observed"), "2");
        assert_eq!(kv(&text, "intended_status"), "all_observed");
        assert_eq!(kv(&text, "intended_missing"), "0");
    }

    #[test]
    fn summary_partial_when_some_intended_missing() {
        let path = scratch("sum-partial").join("summary.txt");
        let result = classify_result(
            2,
            1,
            vec![snp_mut(200, "C", MutationClass::ScatteredSnv)],
            0,
        );
        write_summary(&path, &result, true, "cas9").unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(kv(&text, "intended_declared"), "2");
        assert_eq!(kv(&text, "intended_observed"), "1");
        assert_eq!(kv(&text, "intended_status"), "partial");
        assert_eq!(kv(&text, "intended_missing"), "1");
    }

    #[test]
    fn summary_none_observed_when_declared_positive_and_zero_hits() {
        let path = scratch("sum-none").join("summary.txt");
        let result = classify_result(
            1,
            0,
            vec![snp_mut(200, "C", MutationClass::ScatteredSnv)],
            0,
        );
        write_summary(&path, &result, true, "cas9").unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(kv(&text, "intended_declared"), "1");
        assert_eq!(kv(&text, "intended_observed"), "0");
        assert_eq!(kv(&text, "intended_status"), "none_observed");
        assert_eq!(kv(&text, "intended_missing"), "1");
    }

    #[test]
    fn sub_tsv_writes_expected_coordinates_and_alleles() {
        let refs = vec![RefContig {
            name: "chr".into(),
            seq: b"ACGTACGT".to_vec(),
        }];
        let e = GdEntry::sub(1, "chr", 3, 2, "TT");
        let (s, a, b) = coords(&e);
        assert_eq!(s, "chr");
        assert_eq!(a, 3);
        assert_eq!(b, 4);
        let (ref_al, alt_al) = alleles(&e, &refs);
        assert_eq!(ref_al, "GT");
        assert_eq!(alt_al, "TT");
    }

    #[test]
    fn summary_empty_intended_file_is_na() {
        let path = scratch("sum-empty").join("summary.txt");
        let result = classify_result(0, 0, vec![], 0);
        write_summary(&path, &result, true, "cas9").unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(kv(&text, "intended_provided"), "yes");
        assert_eq!(kv(&text, "intended_declared"), "0");
        assert_eq!(kv(&text, "intended_observed"), "0");
        assert_eq!(kv(&text, "intended_status"), "NA");
        assert_eq!(kv(&text, "intended_missing"), "NA");
    }

    #[test]
    fn publish_failure_restores_the_previous_complete_product_set() {
        let path = publication_scratch("publish-rollback");
        let old = publication_audit("old-sentinel");
        let replacement = publication_audit("new-sentinel");
        write_product_outputs(&old, SchemaVersion::V2, &path).expect("publish old product set");
        let expected = product_bytes(&path);

        let fault = inject_fault(PublicationFault::PublishAt(1));
        let result = write_product_outputs(&replacement, SchemaVersion::V2, &path);
        drop(fault);

        assert!(matches!(result, Err(ProductOutputError::Io { .. })));
        assert_eq!(product_bytes(&path), expected);
        assert!(backup_paths(&path).is_empty());
        let visible = product_bytes(&path);
        assert!(visible
            .values()
            .all(|contents| !String::from_utf8_lossy(contents).contains("new-sentinel")));
    }

    #[test]
    fn restore_failure_retains_the_complete_backup_for_manual_recovery() {
        let path = publication_scratch("restore-rollback");
        let old = publication_audit("old-sentinel");
        let replacement = publication_audit("new-sentinel");
        write_product_outputs(&old, SchemaVersion::V2, &path).expect("publish old product set");
        let expected = product_bytes(&path);

        let fault = inject_fault(PublicationFault::PublishAtThenRestoreAt {
            publish: 1,
            restore: 0,
        });
        let result = write_product_outputs(&replacement, SchemaVersion::V2, &path);
        drop(fault);

        let backup_path = match result {
            Err(ProductOutputError::Rollback { backup_path, .. }) => backup_path,
            other => panic!("expected typed rollback error, got {other:?}"),
        };
        assert!(backup_path.is_dir());
        assert_eq!(product_bytes(&backup_path), expected);
        for name in PRODUCT_OUTPUT_NAMES {
            let path = path.join(name);
            assert!(
                !path.exists()
                    || !String::from_utf8_lossy(&std::fs::read(path).expect("read visible output"))
                        .contains("new-sentinel")
            );
        }
    }
}
