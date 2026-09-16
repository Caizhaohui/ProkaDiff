mod markdown;
mod tables;

pub use markdown::write_markdown_report;
pub use tables::{write_edit_outcomes_tsv, write_post_edit_variants_tsv, write_provenance_tsv};

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use prokadiff_classify::{ClassifiedMutation, ClassifyResult, MutationClass, RefContig};
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
        let events: usize = assessments.iter().map(|a| a.matched_event_ids.len()).sum();
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
    if e.kind == GdKind::Del || e.kind == GdKind::Sub {
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
    use prokadiff_classify::{ClassifiedMutation, MutationClass, RefContig};
    use prokadiff_gd::GdEntry;

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
            intended_declared: declared,
            starter_vs_ref,
            intended_edit_assessments: None,
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
}
