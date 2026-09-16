use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use prokadiff_classify::{AnalysisProvenance, AnnotatedVariant, IntendedEditAssessment, RefContig};
use prokadiff_gd::GdKind;

/// Write the per-edit outcome verification audit table (`edit_outcomes.tsv`).
pub fn write_edit_outcomes_tsv(
    path: impl AsRef<Path>,
    assessments: &[IntendedEditAssessment],
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "edit_id\tkind\tseq_id\texpected_start\texpected_end\tstatus\tmatched_event_ids\tleft_boundary_status\tright_boundary_status\texpected_size\tobserved_size\tunexpected_events\tnotes"
    )?;

    for a in assessments {
        let matched_ids_str = if a.matched_event_ids.is_empty() {
            "NONE".to_string()
        } else {
            a.matched_event_ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };

        let left_status = match &a.left_boundary {
            Some(b) => {
                if b.passed {
                    "PASS".to_string()
                } else {
                    format!("FAIL({}bp)", b.diff_bp)
                }
            }
            None => "NA".to_string(),
        };

        let right_status = match &a.right_boundary {
            Some(b) => {
                if b.passed {
                    "PASS".to_string()
                } else {
                    format!("FAIL({}bp)", b.diff_bp)
                }
            }
            None => "NA".to_string(),
        };

        let exp_size_str = a
            .expected_size
            .map(|s| s.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let obs_size_str = a
            .observed_size
            .map(|s| s.to_string())
            .unwrap_or_else(|| "NA".to_string());

        let unexpected_str = if a.unexpected_event_ids.is_empty() {
            "NONE".to_string()
        } else {
            a.unexpected_event_ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };

        let notes_str = if a.notes.is_empty() {
            "NONE".to_string()
        } else {
            a.notes.join("; ")
        };

        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            a.edit_id,
            a.kind,
            a.seq_id,
            a.expected_start,
            a.expected_end,
            a.status.as_str().to_ascii_uppercase(),
            matched_ids_str,
            left_status,
            right_status,
            exp_size_str,
            obs_size_str,
            unexpected_str,
            notes_str
        )?;
    }

    w.flush()?;
    Ok(())
}

/// Write the unified multi-dimensionally annotated variants table (`post_edit_variants.tsv`).
pub fn write_post_edit_variants_tsv(
    path: impl AsRef<Path>,
    variants: &[AnnotatedVariant],
    refs: &[RefContig],
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "variant_id\tseq_id\tposition\tend\tgd_type\tref\talt\torigin_status\tintended_relation\tsize_class\treview_priority\tguide_relation\tguide_mismatches\tpam\tdist_to_site\tmobile_element\tevidence\tlegacy_class"
    )?;

    for v in variants {
        let (seq_id, pos, end) = variant_coords(v);
        let (ref_allele, alt_allele) = variant_alleles(v, refs);

        let (guide_rel_str, mismatches_str, pam_str, dist_str) = match &v.guide_relation {
            prokadiff_classify::GuideRelation::None => {
                ("NONE", "NA".to_string(), "NA".to_string(), "NA".to_string())
            }
            prokadiff_classify::GuideRelation::OnTarget => (
                "ON_TARGET",
                "0".to_string(),
                "NA".to_string(),
                "0".to_string(),
            ),
            prokadiff_classify::GuideRelation::CandidateOffTarget {
                spacer_mismatches,
                pam,
                distance_to_site,
                ..
            } => (
                "CANDIDATE_OFF_TARGET",
                spacer_mismatches.to_string(),
                pam.clone(),
                distance_to_site.to_string(),
            ),
        };

        let mob_str = match &v.mobile_element_relation {
            Some(m) => {
                let name = m.element_name.as_deref().unwrap_or("IS");
                let tsd = m.target_site_duplication.as_deref().unwrap_or("none");
                format!("{name}(TSD={tsd})")
            }
            None => "NONE".to_string(),
        };

        let legacy_str = v.legacy_class.map(|c| c.as_str()).unwrap_or("none");

        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            v.variant_id,
            seq_id,
            pos,
            end,
            v.entry.kind.as_str(),
            ref_allele,
            alt_allele,
            v.origin_status.as_str(),
            v.intended_relation.as_str(),
            v.size_class.as_str(),
            v.review_priority.as_str(),
            guide_rel_str,
            mismatches_str,
            pam_str,
            dist_str,
            mob_str,
            v.evidence.format_brief(),
            legacy_str
        )?;
    }

    w.flush()?;
    Ok(())
}

/// Write technical provenance information (`provenance.tsv`).
pub fn write_provenance_tsv(
    path: impl AsRef<Path>,
    provenance: &AnalysisProvenance,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "key\tvalue")?;
    writeln!(w, "prokadiff_version\t{}", provenance.prokadiff_version)?;
    writeln!(w, "git_commit\t{}", provenance.git_commit)?;
    writeln!(
        w,
        "reference_sha256\t{}",
        provenance.reference_sha256.as_deref().unwrap_or("NA")
    )?;
    writeln!(
        w,
        "bowtie2_version\t{}",
        provenance.bowtie2_version.as_deref().unwrap_or("NA")
    )?;
    writeln!(
        w,
        "offtarget_search_status\t{}",
        provenance.offtarget_search_status
    )?;
    writeln!(w, "cfd_scoring_status\t{}", provenance.cfd_scoring_status)?;
    writeln!(w, "hsu_scoring_status\t{}", provenance.hsu_scoring_status)?;
    writeln!(w, "bulge_search_status\t{}", provenance.bulge_search_status)?;
    writeln!(w, "run_timestamp\t{}", provenance.run_timestamp)?;

    w.flush()?;
    Ok(())
}

fn variant_coords(v: &AnnotatedVariant) -> (String, u64, u64) {
    let e = &v.entry;
    match e.kind {
        GdKind::Snp | GdKind::Ins => {
            let s = e.seq_id().unwrap_or(".").to_string();
            let p = e.position().unwrap_or(0);
            (s, p, p)
        }
        GdKind::Del | GdKind::Sub => {
            let s = e.fields.first().cloned().unwrap_or_else(|| ".".into());
            let p: u64 = e.fields.get(1).and_then(|x| x.parse().ok()).unwrap_or(0);
            let size: u64 = e.fields.get(2).and_then(|x| x.parse().ok()).unwrap_or(1);
            let end = p.saturating_add(size.saturating_sub(1));
            (s, p, end)
        }
        GdKind::Jc => {
            let s = e.fields.first().cloned().unwrap_or_else(|| ".".into());
            let p: u64 = e.fields.get(1).and_then(|x| x.parse().ok()).unwrap_or(0);
            (s, p, p)
        }
        _ => {
            let s = e.seq_id().unwrap_or(".").to_string();
            let p = e.position().unwrap_or(0);
            (s, p, p)
        }
    }
}

pub(crate) fn variant_alleles(v: &AnnotatedVariant, refs: &[RefContig]) -> (String, String) {
    let e = &v.entry;
    match e.kind {
        GdKind::Snp => {
            let seq_id = e.seq_id().unwrap_or("");
            let pos = e.position().unwrap_or(0);
            let ref_b = if pos > 0 {
                refs.iter()
                    .find(|c| c.name == seq_id)
                    .and_then(|c| c.seq.get((pos - 1) as usize))
                    .map(|&b| (b as char).to_ascii_uppercase().to_string())
                    .unwrap_or_else(|| "N".into())
            } else {
                "N".into()
            };
            let alt = e.fields.get(2).cloned().unwrap_or_else(|| ".".into());
            (ref_b, alt)
        }
        GdKind::Sub => {
            let ref_seq = e.fields.get(2).cloned().unwrap_or_else(|| ".".into());
            let alt_seq = e.fields.get(3).cloned().unwrap_or_else(|| ".".into());
            (ref_seq, alt_seq)
        }
        GdKind::Ins => {
            let alt = e.fields.get(2).cloned().unwrap_or_else(|| ".".into());
            (".".into(), alt)
        }
        GdKind::Del => {
            let size = e.fields.get(2).cloned().unwrap_or_else(|| ".".into());
            (format!("{size}bp"), "none".into())
        }
        GdKind::Mob => {
            let name = e.fields.get(2).cloned().unwrap_or_else(|| "IS".into());
            ("none".into(), name)
        }
        GdKind::Jc => {
            let s2 = e.fields.get(3).cloned().unwrap_or_else(|| ".".into());
            let p2 = e.fields.get(4).cloned().unwrap_or_else(|| ".".into());
            ("none".into(), format!("{s2}:{p2}"))
        }
        _ => (".".into(), ".".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use prokadiff_classify::IntendedEditStatus;

    #[test]
    fn test_write_edit_outcomes_tsv() {
        let temp_file = std::env::temp_dir().join("test_edit_outcomes.tsv");
        let assessments = vec![IntendedEditAssessment {
            edit_id: "edit_1".into(),
            kind: "del".into(),
            seq_id: "chr".into(),
            expected_start: 100,
            expected_end: 200,
            status: IntendedEditStatus::Complete,
            matched_event_ids: vec![12],
            left_boundary: None,
            right_boundary: None,
            expected_size: Some(101),
            observed_size: Some(101),
            unexpected_event_ids: vec![],
            notes: vec!["Both junctions confirmed".into()],
        }];

        write_edit_outcomes_tsv(&temp_file, &assessments).unwrap();
        let content = std::fs::read_to_string(&temp_file).unwrap();
        assert!(content.contains("edit_id\tkind\tseq_id"));
        assert!(content.contains("edit_1\tdel\tchr\t100\t200\tCOMPLETE\t12\tNA\tNA\t101\t101\tNONE\tBoth junctions confirmed"));
        let _ = std::fs::remove_file(temp_file);
    }
}
