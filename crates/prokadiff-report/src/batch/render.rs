//! Multi-clone batch audit report and variants matrix rendering.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use super::recurrence::{BatchVariantRecord, RecurrenceStatus};

/// Summary metrics for an individual clone in the batch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CloneBatchSummary {
    pub clone_id: String,
    pub r1: String,
    pub r2: Option<String>,
    pub intended_declared: usize,
    pub intended_complete: usize,
    pub intended_outcome_display: String,
    pub total_differential_events: usize,
    pub unintended_mutations: usize,
    pub class_1_near_homolog: usize,
    pub class_2_scattered_snv: usize,
    pub class_3_structural: usize,
}

impl CloneBatchSummary {
    pub fn from_audit(
        clone_id: String,
        r1: String,
        r2: Option<String>,
        audit: &prokadiff_classify::AuditResult,
    ) -> Self {
        let declared = audit.summary.intended_declared.unwrap_or(0);
        let complete = audit.summary.intended_edits_complete.unwrap_or(0);
        let intended_outcome_display = match audit.summary.intended_status {
            None => {
                if declared == 0 {
                    "Not specified".to_string()
                } else {
                    format!("None observed (0/{declared})")
                }
            }
            Some(prokadiff_classify::IntendedSummaryStatus::AllObserved) => {
                format!("Complete ({complete}/{declared})")
            }
            Some(prokadiff_classify::IntendedSummaryStatus::NoneObserved) => {
                format!("None observed (0/{declared})")
            }
            Some(prokadiff_classify::IntendedSummaryStatus::Partial) => {
                format!("Partial ({complete}/{declared})")
            }
        };

        Self {
            clone_id,
            r1,
            r2,
            intended_declared: declared,
            intended_complete: complete,
            intended_outcome_display,
            total_differential_events: audit.event_index.len(),
            unintended_mutations: audit.summary.unintended_count,
            class_1_near_homolog: audit.summary.near_homolog_count,
            class_2_scattered_snv: audit.summary.scattered_snv_count,
            class_3_structural: audit.summary.structural_count,
        }
    }
}

/// Extract clone variant observations from an individual clone's AuditResult.
pub fn extract_clone_observations(
    clone_id: &str,
    audit: &prokadiff_classify::AuditResult,
    refs: &[prokadiff_classify::RefContig],
) -> Vec<super::recurrence::CloneVariantObservation> {
    let mut obs = Vec::with_capacity(audit.variants.len());
    for v in &audit.variants {
        let (seq_id, pos, end) = crate::tables::variant_coords(v);
        let (ref_allele, alt_allele) = crate::tables::variant_alleles(v, refs);
        let e = &v.entry;

        let canonical_event = v
            .event_id
            .as_ref()
            .and_then(|id| audit.event_index.get(id))
            .map(|diff| diff.canonical.clone())
            .or_else(|| prokadiff_classify::CanonicalEvent::from_gd_entry(e, refs).ok());

        let (canonical_key, side2_seq_id, side2_position, seq_id_eff, pos_eff, end_eff, ref_eff, alt_eff) =
            if let Some(ref can) = canonical_event {
                let key = can.to_preimage_bytes();
                match can {
                    prokadiff_classify::CanonicalEvent::Jc {
                        side_1,
                        side_2,
                        overlap,
                    } => {
                        let s1_id = side_1.contig.seq_id.clone();
                        let s1_pos = side_1.position;
                        let s2_id = side_2.contig.seq_id.clone();
                        let s2_pos = side_2.position;
                        let alt = format!("{s2_id}:{s2_pos}:overlap={overlap}");
                        (
                            key,
                            Some(s2_id),
                            Some(s2_pos),
                            s1_id,
                            s1_pos,
                            s1_pos,
                            "none".to_string(),
                            alt,
                        )
                    }
                    prokadiff_classify::CanonicalEvent::Con {
                        contig,
                        position,
                        size,
                        region,
                    } => {
                        let cid = contig.seq_id.clone();
                        let p = *position;
                        let end_pos = p.saturating_add(size.saturating_sub(1));
                        (
                            key,
                            None,
                            None,
                            cid,
                            p,
                            end_pos,
                            format!("{size}bp"),
                            region.clone(),
                        )
                    }
                    prokadiff_classify::CanonicalEvent::Mob {
                        contig,
                        position,
                        repeat_name,
                        ..
                    } => {
                        let cid = contig.seq_id.clone();
                        let p = *position;
                        (
                            key,
                            None,
                            None,
                            cid,
                            p,
                            p,
                            "none".to_string(),
                            repeat_name.clone(),
                        )
                    }
                    prokadiff_classify::CanonicalEvent::Amp {
                        contig,
                        position,
                        size,
                        new_copy_number,
                    } => {
                        let cid = contig.seq_id.clone();
                        let p = *position;
                        let end_pos = p.saturating_add(size.saturating_sub(1));
                        (
                            key,
                            None,
                            None,
                            cid,
                            p,
                            end_pos,
                            format!("{size}bp"),
                            format!("{new_copy_number}x"),
                        )
                    }
                    _ => (
                        key,
                        None,
                        None,
                        seq_id,
                        pos,
                        end,
                        ref_allele,
                        alt_allele,
                    ),
                }
            } else {
                let (s2_id, s2_pos) = if e.kind == prokadiff_gd::GdKind::Jc {
                    (
                        e.fields.get(3).cloned(),
                        e.fields.get(4).and_then(|x| x.parse().ok()),
                    )
                } else {
                    (None, None)
                };
                let key = format!("{}:{}:{}:{}:{}:{}", e.kind.as_str(), seq_id, pos, end, ref_allele, alt_allele).into_bytes();
                (key, s2_id, s2_pos, seq_id, pos, end, ref_allele, alt_allele)
            };

        let event_id = match &v.event_id {
            Some(ev) => ev.as_str().to_string(),
            None => v.variant_id.clone(),
        };

        obs.push(super::recurrence::CloneVariantObservation {
            clone_id: clone_id.to_string(),
            event_id,
            seq_id: seq_id_eff,
            position: pos_eff,
            end: end_eff,
            gd_type: v.entry.kind.as_str().to_string(),
            ref_allele: ref_eff,
            alt_allele: alt_eff,
            side2_seq_id,
            side2_position,
            canonical_key,
        });
    }
    obs
}

/// Comprehensive summary of a multi-clone batch audit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BatchAuditSummary {
    pub manifest_path: String,
    pub starter_sample_id: String,
    pub starter_r1: String,
    pub starter_r2: Option<String>,
    pub reference_name: String,
    pub clones: Vec<CloneBatchSummary>,
    pub variants_matrix: Vec<BatchVariantRecord>,
}

/// Render the cross-clone variant recurrence matrix to a writer.
pub fn render_batch_variants_tsv(
    w: &mut impl Write,
    records: &[BatchVariantRecord],
    clone_ids: &[String],
) -> std::io::Result<()> {
    // Header
    write!(
        w,
        "batch_variant_id\tseq_id\tposition\tend\tgd_type\tref\talt\tside2_seq_id\tside2_position\trecurrence\tcount\ttotal"
    )?;
    for cid in clone_ids {
        write!(w, "\tclone_{cid}")?;
    }
    writeln!(w)?;

    let total_clones = clone_ids.len();

    // Data rows
    for r in records {
        let (count, total) = match r.recurrence {
            RecurrenceStatus::CloneSpecific => (1, total_clones),
            RecurrenceStatus::SharedSubset { count, total } => (count, total),
            RecurrenceStatus::SharedAll { total } => (total, total),
        };

        let ref_str = if r.ref_allele.is_empty() {
            "."
        } else {
            &r.ref_allele
        };
        let alt_str = if r.alt_allele.is_empty() {
            "."
        } else {
            &r.alt_allele
        };
        let side2_id = r.side2_seq_id.as_deref().unwrap_or(".");
        let side2_pos = r
            .side2_position
            .map(|p| p.to_string())
            .unwrap_or_else(|| ".".to_string());

        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.batch_variant_id,
            r.seq_id,
            r.position,
            r.end,
            r.gd_type,
            ref_str,
            alt_str,
            side2_id,
            side2_pos,
            r.recurrence.as_str(),
            count,
            total,
        )?;

        for cid in clone_ids {
            let val = match r.clone_presence.get(cid) {
                Some(Some(ev_id)) => ev_id.as_str(),
                _ => ".",
            };
            write!(w, "\t{val}")?;
        }
        writeln!(w)?;
    }

    Ok(())
}

/// Write the cross-clone variant recurrence matrix (`batch_variants.tsv`).
pub fn write_batch_variants_tsv(
    path: impl AsRef<Path>,
    records: &[BatchVariantRecord],
    clone_ids: &[String],
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    render_batch_variants_tsv(&mut w, records, clone_ids)?;
    w.flush()?;
    Ok(())
}

/// Render the human-readable multi-clone batch audit report to a writer.
pub fn render_batch_report_md(
    w: &mut impl Write,
    summary: &BatchAuditSummary,
) -> std::io::Result<()> {
    // 1. Title
    writeln!(w, "# ProkaDiff Multi-Clone Batch Audit Report\n")?;

    // 2. Executive Summary
    writeln!(w, "## 1. Executive Summary\n")?;
    writeln!(w, "- **Starter Strain:** `{}`", summary.starter_sample_id)?;
    writeln!(w, "- **Reference Genome:** `{}`", summary.reference_name)?;
    writeln!(w, "- **Total Evaluated Clones:** {}", summary.clones.len())?;
    writeln!(
        w,
        "- **Total Distinct Post-Edit Variant Loci:** {}\n",
        summary.variants_matrix.len()
    )?;

    // Variant recurrence stats
    let total_vars = summary.variants_matrix.len();
    let mut clone_specific = 0usize;
    let mut shared_subset = 0usize;
    let mut shared_all = 0usize;

    for v in &summary.variants_matrix {
        match v.recurrence {
            RecurrenceStatus::CloneSpecific => clone_specific += 1,
            RecurrenceStatus::SharedSubset { .. } => shared_subset += 1,
            RecurrenceStatus::SharedAll { .. } => shared_all += 1,
        }
    }

    writeln!(w, "### Recurrence Distribution")?;
    if total_vars == 0 {
        writeln!(w, "- No post-edit variants detected across evaluated clones.\n")?;
    } else {
        let pct = |n: usize| (n as f64 / total_vars as f64) * 100.0;
        writeln!(
            w,
            "- **Clone-specific (`clone_specific`):** {} loci ({:.1}%)",
            clone_specific,
            pct(clone_specific)
        )?;
        writeln!(
            w,
            "- **Shared across clone subset (`shared_subset`):** {} loci ({:.1}%)",
            shared_subset,
            pct(shared_subset)
        )?;
        writeln!(
            w,
            "- **Shared across all evaluated clones (`shared_all`):** {} loci ({:.1}%)\n",
            shared_all,
            pct(shared_all)
        )?;
    }

    // 3. Batch Cohort Overview
    writeln!(w, "## 2. Batch Cohort Overview\n")?;
    writeln!(
        w,
        "| Clone ID | Intended Outcome | Differential Events | Unintended Mutations | Class 1 (Near-Homolog) | Class 2 (Scattered SNV) | Class 3 (Structural) |"
    )?;
    writeln!(
        w,
        "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
    )?;

    for c in &summary.clones {
        writeln!(
            w,
            "| `{}` | {} | {} | {} | {} | {} | {} |",
            c.clone_id,
            c.intended_outcome_display,
            c.total_differential_events,
            c.unintended_mutations,
            c.class_1_near_homolog,
            c.class_2_scattered_snv,
            c.class_3_structural,
        )?;
    }
    writeln!(w)?;

    // 4. Cross-Clone Variant Recurrence Analysis
    writeln!(w, "## 3. Cross-Clone Variant Recurrence Analysis\n")?;
    let recurrent_vars: Vec<_> = summary
        .variants_matrix
        .iter()
        .filter(|v| !matches!(v.recurrence, RecurrenceStatus::CloneSpecific))
        .collect();

    if recurrent_vars.is_empty() {
        writeln!(
            w,
            "No recurrent variant loci (shared across multiple clones) were observed.\n"
        )?;
    } else {
        writeln!(
            w,
            "| Variant ID | Seq ID | Position | Type | Ref/Alt | Recurrence Status | Present Clones |"
        )?;
        writeln!(
            w,
            "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
        )?;

        for v in recurrent_vars {
            let ref_str = if v.ref_allele.is_empty() {
                "."
            } else {
                &v.ref_allele
            };
            let alt_str = if v.alt_allele.is_empty() {
                "."
            } else {
                &v.alt_allele
            };
            let ref_alt = format!("{ref_str} -> {alt_str}");

            let mut present_clones = Vec::new();
            for c in &summary.clones {
                if let Some(Some(_)) = v.clone_presence.get(&c.clone_id) {
                    present_clones.push(c.clone_id.as_str());
                }
            }
            let clones_str = present_clones.join(", ");

            writeln!(
                w,
                "| `{}` | `{}` | {} | {} | {} | `{}` | {} |",
                v.batch_variant_id,
                v.seq_id,
                v.position,
                v.gd_type,
                ref_alt,
                v.recurrence.as_str(),
                clones_str,
            )?;
        }
        writeln!(w)?;
    }

    // 5. Non-Prescriptive Analysis Limitations & Cautions (Taskbook §66)
    writeln!(
        w,
        "## 4. Non-Prescriptive Analysis Limitations & Cautions\n"
    )?;
    writeln!(
        w,
        "- **Strict Non-Prescriptive Principle (Taskbook §66):** This report provides observational differential audit metrics and variant recurrence tallies across evaluated clones. In accordance with non-prescriptive audit principles, ProkaDiff does not perform clone selection, ranking, or preference evaluation. Downstream clone selection depends on experimental phenotypes, expression titers, and engineering objectives decided by the investigator."
    )?;
    writeln!(
        w,
        "- **Spatial Association vs. Biological Causality:** Spatial association is consistent with a possible guide-dependent event but does not establish cleavage causality. Differential variants identified in single or multiple clones represent genomic observations relative to the sequenced starter strain."
    )?;
    writeln!(
        w,
        "- **Short-Read Structural Resolution:** Variant and junction calling relies on short Illumina read mappings. Extremely repetitive genomic regions or complex structural rearrangements may not be fully resolved."
    )?;
    writeln!(w)?;

    // 6. Technical Provenance
    writeln!(w, "## 5. Technical Provenance\n")?;
    writeln!(w, "| Parameter | Value |")?;
    writeln!(w, "| :--- | :--- |")?;
    writeln!(w, "| Manifest Path | `{}` |", summary.manifest_path)?;
    writeln!(w, "| Starter R1 | `{}` |", summary.starter_r1)?;
    if let Some(r2) = &summary.starter_r2 {
        writeln!(w, "| Starter R2 | `{}` |", r2)?;
    }
    writeln!(
        w,
        "| ProkaDiff Mode | Multi-Clone Engineering Batch Audit Mode (M10) |"
    )?;
    writeln!(w)?;

    Ok(())
}

/// Write the human-readable multi-clone batch audit report (`batch_report.md`).
pub fn write_batch_report_md(
    path: impl AsRef<Path>,
    summary: &BatchAuditSummary,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    render_batch_report_md(&mut w, summary)?;
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_render_batch_variants_tsv_and_report() {
        let clone_ids = vec!["C01".to_string(), "C02".to_string(), "C03".to_string()];

        let mut cp1 = HashMap::new();
        cp1.insert("C01".to_string(), Some("DV1_001".to_string()));
        cp1.insert("C02".to_string(), Some("DV1_002".to_string()));
        cp1.insert("C03".to_string(), Some("DV1_003".to_string()));

        let mut cp2 = HashMap::new();
        cp2.insert("C01".to_string(), Some("DV1_004".to_string()));
        cp2.insert("C02".to_string(), None);
        cp2.insert("C03".to_string(), None);

        let records = vec![
            BatchVariantRecord {
                batch_variant_id: "BV_0001".to_string(),
                seq_id: "chr1".to_string(),
                position: 100,
                end: 100,
                gd_type: "SNP".to_string(),
                ref_allele: "A".to_string(),
                alt_allele: "G".to_string(),
                side2_seq_id: None,
                side2_position: None,
                recurrence: RecurrenceStatus::SharedAll { total: 3 },
                clone_presence: cp1,
            },
            BatchVariantRecord {
                batch_variant_id: "BV_0002".to_string(),
                seq_id: "chr1".to_string(),
                position: 250,
                end: 250,
                gd_type: "DEL".to_string(),
                ref_allele: "T".to_string(),
                alt_allele: "".to_string(),
                side2_seq_id: None,
                side2_position: None,
                recurrence: RecurrenceStatus::CloneSpecific,
                clone_presence: cp2,
            },
        ];

        // 1. Test TSV rendering
        let mut tsv_buf = Vec::new();
        render_batch_variants_tsv(&mut tsv_buf, &records, &clone_ids).unwrap();
        let tsv_content = String::from_utf8(tsv_buf).unwrap();

        assert!(tsv_content.contains("batch_variant_id\tseq_id"));
        assert!(tsv_content.contains("clone_C01\tclone_C02\tclone_C03"));
        assert!(tsv_content.contains("BV_0001\tchr1\t100\t100\tSNP\tA\tG\t.\t.\tshared_all\t3\t3\tDV1_001\tDV1_002\tDV1_003"));
        assert!(tsv_content.contains("BV_0002\tchr1\t250\t250\tDEL\tT\t.\t.\t.\tclone_specific\t1\t3\tDV1_004\t.\t."));

        // 2. Test Report rendering
        let summary = BatchAuditSummary {
            manifest_path: "batch_manifest.tsv".to_string(),
            starter_sample_id: "WT_Starter".to_string(),
            starter_r1: "wt_R1.fq".to_string(),
            starter_r2: None,
            reference_name: "NC_000913.3".to_string(),
            clones: vec![
                CloneBatchSummary {
                    clone_id: "C01".to_string(),
                    r1: "c1_R1.fq".to_string(),
                    r2: None,
                    intended_declared: 1,
                    intended_complete: 1,
                    intended_outcome_display: "Complete (1/1)".to_string(),
                    total_differential_events: 2,
                    unintended_mutations: 1,
                    class_1_near_homolog: 0,
                    class_2_scattered_snv: 1,
                    class_3_structural: 0,
                },
                CloneBatchSummary {
                    clone_id: "C02".to_string(),
                    r1: "c2_R1.fq".to_string(),
                    r2: None,
                    intended_declared: 1,
                    intended_complete: 0,
                    intended_outcome_display: "None observed (0/1)".to_string(),
                    total_differential_events: 1,
                    unintended_mutations: 0,
                    class_1_near_homolog: 0,
                    class_2_scattered_snv: 0,
                    class_3_structural: 0,
                },
                CloneBatchSummary {
                    clone_id: "C03".to_string(),
                    r1: "c3_R1.fq".to_string(),
                    r2: None,
                    intended_declared: 0,
                    intended_complete: 0,
                    intended_outcome_display: "Not specified".to_string(),
                    total_differential_events: 1,
                    unintended_mutations: 0,
                    class_1_near_homolog: 0,
                    class_2_scattered_snv: 0,
                    class_3_structural: 0,
                },
            ],
            variants_matrix: records,
        };

        let mut report_buf = Vec::new();
        render_batch_report_md(&mut report_buf, &summary).unwrap();
        let report_content = String::from_utf8(report_buf).unwrap();

        assert!(report_content.contains("# ProkaDiff Multi-Clone Batch Audit Report"));
        assert!(report_content.contains("Strict Non-Prescriptive Principle (Taskbook §66)"));
        assert!(report_content.contains("Spatial association is consistent with a possible guide-dependent event but does not establish cleavage causality."));
        assert!(report_content.contains("None observed (0/1)"));
        assert!(report_content.contains("BV_0001"));
        assert!(!report_content.contains("best"));
        assert!(!report_content.contains("safest"));
    }

    #[test]
    fn test_extract_clone_observations_canonical_jc_swap_and_con() {
        use prokadiff_classify::{
            AnnotatedVariant, AuditResult, AuditSummary, CandidateSearchStatus, GuideRelation,
            IntendedRelation, OriginStatus, RefContig, ReviewPriority, SampleMetadata, SizeClass,
        };
        use prokadiff_gd::{GdEntry, GdKind};

        let refs = vec![
            RefContig {
                name: "chr1".to_string(),
                seq: vec![b'A'; 2000],
            },
            RefContig {
                name: "chr2".to_string(),
                seq: vec![b'C'; 3000],
            },
        ];

        let mock_meta = SampleMetadata {
            starter_names: vec!["s.fq".into()],
            edited_names: vec!["e.fq".into()],
            reference_names: vec!["ref.fa".into()],
            editor: "dsb".into(),
            spacer: None,
            pam: None,
            threads: 1,
        };

        let mock_provenance = prokadiff_classify::AnalysisProvenance {
            prokadiff_version: "0.1.0".into(),
            git_commit: "test".into(),
            reference_sha256: None,
            bowtie2_version: None,
            offtarget_search_status: "NOT_REQUESTED".into(),
            cfd_scoring_status: "DISABLED".into(),
            hsu_scoring_status: "DISABLED".into(),
            bulge_search_status: "EXACT".into(),
            run_timestamp: "2026-11-04T00:00:00Z".into(),
        };

        let empty_summary = AuditSummary {
            intended_provided: false,
            intended_declared: None,
            intended_edits_complete: None,
            intended_edits_partial: None,
            intended_edits_missing: None,
            intended_edits_unexpected: None,
            intended_events_observed: None,
            intended_status: None,
            intended_missing: None,
            starter_vs_reference: 0,
            unintended_count: 0,
            post_edit_variant_count: 0,
            structural_count: 0,
            near_homolog_count: 0,
            scattered_snv_count: 0,
            mobile_element_count: 0,
            non_mobile_structural_count: 0,
            candidate_offtarget_count: 0,
            unassociated_small_count: 0,
            high_attention_count: 0,
            review_count: 0,
            info_count: 0,
            association_count: 0,
            candidate_search_status: CandidateSearchStatus::NotPerformed,
        };

        // Clone 1 has JC: chr1 100 1 chr2 200 -1 0, and CON: chr1 500 20 chr1:1000-1019
        let entry_jc1 = GdEntry {
            kind: GdKind::Jc,
            id: 1,
            parent_ids: Vec::new(),
            fields: vec![
                "chr1".into(),
                "100".into(),
                "1".into(),
                "chr2".into(),
                "200".into(),
                "-1".into(),
                "0".into(),
            ],
            attrs: std::collections::BTreeMap::new(),
        };

        let entry_con1 = GdEntry {
            kind: GdKind::Con,
            id: 2,
            parent_ids: Vec::new(),
            fields: vec!["chr1".into(), "500".into(), "20".into(), "chr1:1000-1019".into()],
            attrs: std::collections::BTreeMap::new(),
        };

        let audit1 = AuditResult {
            sample: mock_meta.clone(),
            intended_provided: false,
            intended_edits: vec![],
            intended_observed_event_ids: vec![],
            starter_vs_reference: 0,
            hypothesis_enabled: false,
            event_index: std::collections::BTreeMap::new(),
            event_display: std::collections::BTreeMap::new(),
            unintended: vec![],
            variants: vec![
                AnnotatedVariant {
                    variant_id: "VAR_0001".into(),
                    event_id: None,
                    entry: entry_jc1,
                    origin_status: OriginStatus::PostEditDifferential,
                    size_class: SizeClass::Structural,
                    intended_relation: IntendedRelation::None,
                    guide_relation: GuideRelation::None,
                    mobile_element_relation: None,
                    repeat_relation: None,
                    gene_annotation: None,
                    evidence: prokadiff_classify::EvidenceSummary::unknown(),
                    review_priority: ReviewPriority::Info,
                    legacy_class: None,
                    hypothesis: None,
                },
                AnnotatedVariant {
                    variant_id: "VAR_0002".into(),
                    event_id: None,
                    entry: entry_con1,
                    origin_status: OriginStatus::PostEditDifferential,
                    size_class: SizeClass::Structural,
                    intended_relation: IntendedRelation::None,
                    guide_relation: GuideRelation::None,
                    mobile_element_relation: None,
                    repeat_relation: None,
                    gene_annotation: None,
                    evidence: prokadiff_classify::EvidenceSummary::unknown(),
                    review_priority: ReviewPriority::Info,
                    legacy_class: None,
                    hypothesis: None,
                },
            ],
            guide_sites: vec![],
            variant_site_links: vec![],
            associations: vec![],
            candidate_search_status: CandidateSearchStatus::NotPerformed,
            summary: empty_summary.clone(),
            provenance: mock_provenance.clone(),
            subclonal_enabled: false,
            subclonal_variants: vec![],
        };

        // Clone 2 has endpoint-swapped JC: chr2 200 -1 chr1 100 1 0, and CON: chr1 500 30 chr2:2000-2029 (diff size & region)
        let entry_jc2 = GdEntry {
            kind: GdKind::Jc,
            id: 1,
            parent_ids: Vec::new(),
            fields: vec![
                "chr2".into(),
                "200".into(),
                "-1".into(),
                "chr1".into(),
                "100".into(),
                "1".into(),
                "0".into(),
            ],
            attrs: std::collections::BTreeMap::new(),
        };

        let entry_con2 = GdEntry {
            kind: GdKind::Con,
            id: 2,
            parent_ids: Vec::new(),
            fields: vec!["chr1".into(), "500".into(), "30".into(), "chr2:2000-2029".into()],
            attrs: std::collections::BTreeMap::new(),
        };

        let audit2 = AuditResult {
            sample: mock_meta,
            intended_provided: false,
            intended_edits: vec![],
            intended_observed_event_ids: vec![],
            starter_vs_reference: 0,
            hypothesis_enabled: false,
            event_index: std::collections::BTreeMap::new(),
            event_display: std::collections::BTreeMap::new(),
            unintended: vec![],
            variants: vec![
                AnnotatedVariant {
                    variant_id: "VAR_0001".into(),
                    event_id: None,
                    entry: entry_jc2,
                    origin_status: OriginStatus::PostEditDifferential,
                    size_class: SizeClass::Structural,
                    intended_relation: IntendedRelation::None,
                    guide_relation: GuideRelation::None,
                    mobile_element_relation: None,
                    repeat_relation: None,
                    gene_annotation: None,
                    evidence: prokadiff_classify::EvidenceSummary::unknown(),
                    review_priority: ReviewPriority::Info,
                    legacy_class: None,
                    hypothesis: None,
                },
                AnnotatedVariant {
                    variant_id: "VAR_0002".into(),
                    event_id: None,
                    entry: entry_con2,
                    origin_status: OriginStatus::PostEditDifferential,
                    size_class: SizeClass::Structural,
                    intended_relation: IntendedRelation::None,
                    guide_relation: GuideRelation::None,
                    mobile_element_relation: None,
                    repeat_relation: None,
                    gene_annotation: None,
                    evidence: prokadiff_classify::EvidenceSummary::unknown(),
                    review_priority: ReviewPriority::Info,
                    legacy_class: None,
                    hypothesis: None,
                },
            ],
            guide_sites: vec![],
            variant_site_links: vec![],
            associations: vec![],
            candidate_search_status: CandidateSearchStatus::NotPerformed,
            summary: empty_summary,
            provenance: mock_provenance,
            subclonal_enabled: false,
            subclonal_variants: vec![],
        };

        // Extract observations using production extractor
        let obs1 = extract_clone_observations("C01", &audit1, &refs);
        let obs2 = extract_clone_observations("C02", &audit2, &refs);

        // 1. Verify JC canonicalization: both produce identical canonical_key and canonical coordinates
        assert_eq!(obs1[0].canonical_key, obs2[0].canonical_key);
        assert_eq!(obs1[0].seq_id, "chr1");
        assert_eq!(obs1[0].position, 100);
        assert_eq!(obs1[0].side2_seq_id, Some("chr2".into()));
        assert_eq!(obs1[0].side2_position, Some(200));
        assert_eq!(obs2[0].seq_id, "chr1");
        assert_eq!(obs2[0].position, 100);
        assert_eq!(obs2[0].side2_seq_id, Some("chr2".into()));
        assert_eq!(obs2[0].side2_position, Some(200));

        // 2. Verify CON distinction: different sizes and donor regions produce different canonical_keys
        assert_ne!(obs1[1].canonical_key, obs2[1].canonical_key);

        // 3. Align batch variants across both clones
        let mut all_obs = obs1;
        all_obs.extend(obs2);
        let aligned = crate::batch::align_batch_variants(&["C01".into(), "C02".into()], &all_obs);

        // Expect: 1 merged JC (SharedAll), and 2 separate CONs (CloneSpecific)
        assert_eq!(aligned.len(), 3);
        let jc = aligned.iter().find(|r| r.gd_type == "JC").unwrap();
        assert_eq!(jc.recurrence, RecurrenceStatus::SharedAll { total: 2 });
        assert_eq!(jc.clone_presence.get("C01").unwrap(), &Some("VAR_0001".into()));
        assert_eq!(jc.clone_presence.get("C02").unwrap(), &Some("VAR_0001".into()));

        let cons: Vec<_> = aligned.iter().filter(|r| r.gd_type == "CON").collect();
        assert_eq!(cons.len(), 2);
        assert_eq!(cons[0].recurrence, RecurrenceStatus::CloneSpecific);
        assert_eq!(cons[1].recurrence, RecurrenceStatus::CloneSpecific);
    }
}
