use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use prokadiff_classify::{
    AnnotatedVariant, AuditResult, GuideRelation, IntendedEditStatus, OriginStatus, RefContig,
    ReviewPriority, SizeClass,
};
use prokadiff_gd::GdKind;

/// Render the complete human-readable Genome Audit Report (`report.md`).
pub fn write_markdown_report(
    path: impl AsRef<Path>,
    audit: &AuditResult,
    refs: &[RefContig],
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);

    // 1. Report Title
    writeln!(w, "# ProkaDiff Genome Audit Report\n")?;

    // 2. Executive Summary
    write_executive_summary(&mut w, audit)?;

    // 3. Sample Context
    write_sample_context(&mut w, audit)?;

    // 4. Intended Edit Assessment
    write_intended_assessment(&mut w, audit)?;

    // 5. Genome-wide Findings Overview
    write_differential_summary(&mut w, audit)?;

    // 6. Structural & Mobile Element Events
    write_structural_events(&mut w, audit, refs)?;

    // 7. Candidate Guide-dependent Off-target Events
    write_candidate_offtarget_events(&mut w, audit, refs)?;

    // 8. Distal Small Variants
    write_distal_small_variants(&mut w, audit, refs)?;

    // 9. Analysis Limitations
    write_limitations(&mut w, audit)?;

    // 10. Technical Provenance
    write_provenance(&mut w, audit)?;

    w.flush()?;
    Ok(())
}

fn write_executive_summary(w: &mut impl Write, audit: &AuditResult) -> std::io::Result<()> {
    writeln!(w, "## 1. Executive Summary\n")?;

    // Intended edits status line
    if audit.intended_edits.is_empty() {
        writeln!(
            w,
            "**Intended Edit Outcome:** *Not Specified* (Analysis run in untargeted differential mode).\n"
        )?;
    } else {
        let complete = audit
            .intended_edits
            .iter()
            .filter(|a| a.status == IntendedEditStatus::Complete)
            .count();
        let partial = audit
            .intended_edits
            .iter()
            .filter(|a| a.status == IntendedEditStatus::Partial)
            .count();
        let missing = audit
            .intended_edits
            .iter()
            .filter(|a| a.status == IntendedEditStatus::Missing)
            .count();
        let unexpected = audit
            .intended_edits
            .iter()
            .filter(|a| a.status == IntendedEditStatus::UnexpectedStructure)
            .count();

        if complete == audit.intended_edits.len() {
            writeln!(
                w,
                "**Intended Edit Outcome: PASS** — All declared intended edits ({}) were confirmed with expected genomic structures.\n",
                complete
            )?;
        } else if unexpected > 0 {
            writeln!(
                w,
                "**Intended Edit Outcome: ALERT** — {} intended edit(s) exhibited unexpected/aberrant genomic structures at on-target loci.\n",
                unexpected
            )?;
        } else if partial > 0 {
            writeln!(
                w,
                "**Intended Edit Outcome: PARTIAL** — {} edit(s) complete, {} edit(s) partial, {} edit(s) missing.\n",
                complete, partial, missing
            )?;
        } else {
            writeln!(
                w,
                "**Intended Edit Outcome: MISSING** — 0 intended edits were confirmed across declared loci.\n"
            )?;
        }
    }

    // Post-edit differential variants count
    let post_edit_vars: Vec<&AnnotatedVariant> = audit
        .variants
        .iter()
        .filter(|v| v.origin_status == OriginStatus::PostEditDifferential)
        .collect();

    let n_total = post_edit_vars.len();
    let n_mob = post_edit_vars
        .iter()
        .filter(|v| v.entry.kind == GdKind::Mob)
        .count();
    let n_struct = post_edit_vars
        .iter()
        .filter(|v| v.size_class == SizeClass::Structural && v.entry.kind != GdKind::Mob)
        .count();
    let n_candidate_ot = post_edit_vars
        .iter()
        .filter(|v| matches!(v.guide_relation, GuideRelation::CandidateOffTarget { .. }))
        .count();
    let n_small = post_edit_vars
        .iter()
        .filter(|v| {
            v.size_class == SizeClass::Small && matches!(v.guide_relation, GuideRelation::None)
        })
        .count();

    writeln!(w, "### Key Findings Breakdown")?;
    writeln!(
        w,
        "- **Total post-edit differential variants:** {}",
        n_total
    )?;
    writeln!(w, "- **Mobile element (IS) insertions:** {}", n_mob)?;
    writeln!(
        w,
        "- **Genomic rearrangements & large structural deletions:** {}",
        n_struct
    )?;
    writeln!(
        w,
        "- **Candidate guide-dependent off-target events:** {}",
        n_candidate_ot
    )?;
    writeln!(w, "- **Distal / collateral small variants:** {}\n", n_small)?;

    if n_total == 0 {
        writeln!(
            w,
            "> **Audit Conclusion:** No additional post-edit differential variants were detected outside declared edits.\n"
        )?;
    } else {
        let n_high = post_edit_vars
            .iter()
            .filter(|v| v.review_priority == ReviewPriority::HighAttention)
            .count();
        let n_rev = post_edit_vars
            .iter()
            .filter(|v| v.review_priority == ReviewPriority::Review)
            .count();

        writeln!(
            w,
            "> **Audit Review Recommendation:** {} high-attention structural or near-target event(s) and {} review-priority variant(s) require manual inspection.",
            n_high, n_rev
        )?;
        writeln!(
            w,
            "> *ProkaDiff provides observational audit priorities, not biological safety claims. Final strain suitability depends on specific experimental requirements.*\n"
        )?;
    }

    Ok(())
}

fn write_sample_context(w: &mut impl Write, audit: &AuditResult) -> std::io::Result<()> {
    writeln!(w, "## 2. Sample & Experiment Context\n")?;
    writeln!(w, "| Parameter | Declared Value |")?;
    writeln!(w, "| :--- | :--- |")?;
    writeln!(
        w,
        "| Starter strain WGS | `{}` |",
        audit.sample.starter_names.join(", ")
    )?;
    writeln!(
        w,
        "| Edited strain WGS | `{}` |",
        audit.sample.edited_names.join(", ")
    )?;
    writeln!(
        w,
        "| Reference genome | `{}` |",
        audit.sample.reference_names.join(", ")
    )?;
    writeln!(w, "| Editor | `{}` |", audit.sample.editor)?;
    if let Some(sp) = &audit.sample.spacer {
        writeln!(w, "| Guide spacer | `{}` |", sp)?;
    }
    if let Some(pam) = &audit.sample.pam {
        writeln!(w, "| PAM profile | `{}` |", pam)?;
    }
    writeln!(w, "| Parallel threads | `{}` |\n", audit.sample.threads)?;
    Ok(())
}

fn write_intended_assessment(w: &mut impl Write, audit: &AuditResult) -> std::io::Result<()> {
    writeln!(w, "## 3. Intended Edit Assessment\n")?;
    if audit.intended_edits.is_empty() {
        writeln!(
            w,
            "*No intended edits were declared. To verify targeting outcomes, supply an `--intended` TSV table.*\n"
        )?;
        return Ok(());
    }

    writeln!(
        w,
        "| Edit ID | Kind | Locus | Status | Left Boundary | Right Boundary | Exp / Obs Size | Diagnostics |"
    )?;
    writeln!(
        w,
        "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- |"
    )?;

    for a in &audit.intended_edits {
        let left_str = match &a.left_boundary {
            Some(b) if b.passed => "PASS".to_string(),
            Some(b) => format!("FAIL({}bp)", b.diff_bp),
            None => "NA".to_string(),
        };
        let right_str = match &a.right_boundary {
            Some(b) if b.passed => "PASS".to_string(),
            Some(b) => format!("FAIL({}bp)", b.diff_bp),
            None => "NA".to_string(),
        };

        let size_str = format!(
            "{} / {}",
            a.expected_size
                .map(|s| s.to_string())
                .unwrap_or_else(|| "NA".into()),
            a.observed_size
                .map(|s| s.to_string())
                .unwrap_or_else(|| "NA".into())
        );

        let notes = if a.notes.is_empty() {
            "none".to_string()
        } else {
            a.notes.join("; ")
        };

        writeln!(
            w,
            "| **{}** | `{}` | {}:{}-{} | **{}** | {} | {} | {} | {} |",
            a.edit_id,
            a.kind,
            a.seq_id,
            a.expected_start,
            a.expected_end,
            a.status.as_str().to_ascii_uppercase(),
            left_str,
            right_str,
            size_str,
            notes
        )?;
    }
    writeln!(w)?;
    Ok(())
}

fn write_differential_summary(w: &mut impl Write, audit: &AuditResult) -> std::io::Result<()> {
    writeln!(w, "## 4. Genome-wide Differential Summary\n")?;

    let post_edit: Vec<&AnnotatedVariant> = audit
        .variants
        .iter()
        .filter(|v| v.origin_status == OriginStatus::PostEditDifferential)
        .collect();

    if post_edit.is_empty() {
        writeln!(
            w,
            "No post-edit differential variants were detected in the edited strain relative to the starter strain.\n"
        )?;
        return Ok(());
    }

    let mut high_count = 0;
    let mut rev_count = 0;
    let mut info_count = 0;

    for v in &post_edit {
        match v.review_priority {
            ReviewPriority::HighAttention => high_count += 1,
            ReviewPriority::Review => rev_count += 1,
            ReviewPriority::Info => info_count += 1,
        }
    }

    writeln!(
        w,
        "| Priority Category | Count | Primary Causes / Characteristics |"
    )?;
    writeln!(w, "| :--- | :---: | :--- |")?;
    writeln!(
        w,
        "| **HIGH_ATTENTION** | **{}** | Large structural rearrangements, IS transposon insertions, or aberrant target structures |",
        high_count
    )?;
    writeln!(
        w,
        "| **REVIEW** | **{}** | Distal point mutations, small indels, or non-target coding variants |",
        rev_count
    )?;
    writeln!(
        w,
        "| **INFO** | **{}** | Confirmed expected edits and neutral background changes |",
        info_count
    )?;
    writeln!(w)?;
    Ok(())
}

fn write_structural_events(
    w: &mut impl Write,
    audit: &AuditResult,
    _refs: &[RefContig],
) -> std::io::Result<()> {
    writeln!(w, "## 5. Structural and Mobile-element Events\n")?;

    let struct_vars: Vec<&AnnotatedVariant> = audit
        .variants
        .iter()
        .filter(|v| {
            v.origin_status == OriginStatus::PostEditDifferential
                && v.size_class == SizeClass::Structural
        })
        .collect();

    if struct_vars.is_empty() {
        writeln!(
            w,
            "No post-edit differential structural or mobile-element events were detected.\n"
        )?;
        return Ok(());
    }

    writeln!(
        w,
        "| Variant ID | Type | Locus | Description | Evidence | Priority |"
    )?;
    writeln!(w, "| :--- | :---: | :--- | :--- | :---: | :---: |")?;

    for v in struct_vars {
        let locus = format!(
            "{}:{}",
            v.entry.seq_id().unwrap_or("."),
            v.entry.position().unwrap_or(0)
        );

        let desc = match &v.mobile_element_relation {
            Some(m) => {
                let el = m.element_name.as_deref().unwrap_or("IS");
                let fam_part = match &m.family {
                    Some(f)
                        if f != el
                            && !f.to_ascii_uppercase().starts_with(&el.to_ascii_uppercase()) =>
                    {
                        format!(" [{f}]")
                    }
                    _ => String::new(),
                };
                format!(
                    "Mobile element {}{} insertion (TSD={})",
                    el,
                    fam_part,
                    m.target_site_duplication.as_deref().unwrap_or("none")
                )
            }
            None => match v.entry.kind {
                GdKind::Del => {
                    let sz = v.entry.fields.get(2).cloned().unwrap_or_else(|| ".".into());
                    format!("Structural deletion (size: {sz} bp)")
                }
                GdKind::Jc => {
                    let s2 = v.entry.fields.get(3).cloned().unwrap_or_else(|| ".".into());
                    let p2 = v.entry.fields.get(4).cloned().unwrap_or_else(|| ".".into());
                    format!("Novel sequence junction linked to {s2}:{p2}")
                }
                _ => v.entry.kind.as_str().to_string(),
            },
        };

        writeln!(
            w,
            "| **{}** | `{}` | {} | {} | `{}` | **{}** |",
            v.variant_id,
            v.entry.kind.as_str(),
            locus,
            desc,
            v.evidence.format_brief(),
            v.review_priority.as_str()
        )?;
    }
    writeln!(w)?;
    Ok(())
}

fn write_candidate_offtarget_events(
    w: &mut impl Write,
    audit: &AuditResult,
    _refs: &[RefContig],
) -> std::io::Result<()> {
    writeln!(w, "## 6. Candidate Guide-dependent Off-target Events\n")?;

    let ot_vars: Vec<&AnnotatedVariant> = audit
        .variants
        .iter()
        .filter(|v| {
            v.origin_status == OriginStatus::PostEditDifferential
                && matches!(v.guide_relation, GuideRelation::CandidateOffTarget { .. })
        })
        .collect();

    if ot_vars.is_empty() {
        writeln!(
            w,
            "No candidate guide-dependent off-target events were identified within the configured search window.\n"
        )?;
        return Ok(());
    }

    writeln!(
        w,
        "| Variant ID | Locus | Change | Mismatches | PAM | Distance to Site | Review Priority |"
    )?;
    writeln!(w, "| :--- | :--- | :---: | :---: | :---: | :---: | :---: |")?;

    for v in ot_vars {
        let locus = format!(
            "{}:{}",
            v.entry.seq_id().unwrap_or("."),
            v.entry.position().unwrap_or(0)
        );

        let (mismatches, pam, dist) = match &v.guide_relation {
            GuideRelation::CandidateOffTarget {
                spacer_mismatches,
                pam,
                distance_to_site,
                ..
            } => (*spacer_mismatches, pam.as_str(), *distance_to_site),
            _ => (0, "NA", 0),
        };

        let change = match v.entry.kind {
            GdKind::Snp => v
                .entry
                .fields
                .get(2)
                .cloned()
                .unwrap_or_else(|| "SNP".into()),
            _ => v.entry.kind.as_str().to_string(),
        };

        writeln!(
            w,
            "| **{}** | {} | `{}` | {} | `{}` | {} bp | **{}** |",
            v.variant_id,
            locus,
            change,
            mismatches,
            pam,
            dist,
            v.review_priority.as_str()
        )?;
    }

    writeln!(
        w,
        "\n> **Scientific Caution:** Spatial proximity and sequence homology to predicted guide sites indicate candidate association, but do not establish nuclease cleavage causality without experimental controls.\n"
    )?;
    Ok(())
}

fn write_distal_small_variants(
    w: &mut impl Write,
    audit: &AuditResult,
    refs: &[RefContig],
) -> std::io::Result<()> {
    writeln!(w, "## 7. Distal Small Variants\n")?;

    let distal_vars: Vec<&AnnotatedVariant> = audit
        .variants
        .iter()
        .filter(|v| {
            v.origin_status == OriginStatus::PostEditDifferential
                && v.size_class == SizeClass::Small
                && matches!(v.guide_relation, GuideRelation::None)
        })
        .collect();

    if distal_vars.is_empty() {
        writeln!(w, "No distal small variants were detected.\n")?;
        return Ok(());
    }

    writeln!(
        w,
        "{} small post-edit differential variant(s) were detected without nearby predicted guide-homologous sites.\n",
        distal_vars.len()
    )?;

    writeln!(
        w,
        "| Variant ID | Type | Locus | Gene / Region | Reference | Alternate | Evidence | Priority |"
    )?;
    writeln!(
        w,
        "| :--- | :---: | :--- | :--- | :---: | :---: | :---: | :---: |"
    )?;

    for v in distal_vars {
        let locus = format!(
            "{}:{}",
            v.entry.seq_id().unwrap_or("."),
            v.entry.position().unwrap_or(0)
        );

        let gene_str = match &v.gene_annotation {
            Some(g) => {
                let name = g
                    .gene_name
                    .as_deref()
                    .or(g.locus_tag.as_deref())
                    .unwrap_or("unnamed");
                if g.feature_type == "intergenic" {
                    g.consequence
                        .clone()
                        .unwrap_or_else(|| format!("intergenic ({name})"))
                } else {
                    let locus_part = g
                        .locus_tag
                        .as_deref()
                        .map(|l| format!(" ({l})"))
                        .unwrap_or_default();
                    format!("{name}{locus_part} [{}]", g.feature_type)
                }
            }
            None => "-".to_string(),
        };

        let (ref_b, alt_b) = crate::tables::variant_alleles(v, refs);

        writeln!(
            w,
            "| **{}** | `{}` | {} | {} | `{}` | `{}` | `{}` | {} |",
            v.variant_id,
            v.entry.kind.as_str(),
            locus,
            gene_str,
            ref_b,
            alt_b,
            v.evidence.format_brief(),
            v.review_priority.as_str()
        )?;
    }

    writeln!(
        w,
        "\n> **Interpretation:** Distal variants lack sequence and spatial association with the declared editing system. Potential origins include spontaneous culture drift, electroporation/transformation stress, or uncharacterized cellular burden.\n"
    )?;
    Ok(())
}

fn write_limitations(w: &mut impl Write, audit: &AuditResult) -> std::io::Result<()> {
    writeln!(w, "## 8. Analysis Limitations & Cautions\n")?;
    writeln!(
        w,
        "- **Short-Read Structural Resolution:** Variant and junction calling relies on short Illumina read mappings. Extremely repetitive genomic regions or large balanced inversions may not be fully resolved."
    )?;
    writeln!(
        w,
        "- **Observational Association vs. Biological Causality:** Post-edit differential variants cannot be automatically attributed as direct enzymatic consequences of editing without independent controls or multiple clones."
    )?;
    if audit.provenance.cfd_scoring_status.contains("DISABLED")
        || audit.provenance.cfd_scoring_status.contains("UNVALIDATED")
    {
        writeln!(
            w,
            "- **CFD Scoring Unavailable:** CFD off-target cleavage scores were omitted because external-oracle parity validation was not enabled for this analysis."
        )?;
    }
    writeln!(w)?;
    Ok(())
}

fn write_provenance(w: &mut impl Write, audit: &AuditResult) -> std::io::Result<()> {
    writeln!(w, "## 9. Technical Provenance\n")?;
    writeln!(w, "| Parameter | Provenance Detail |")?;
    writeln!(w, "| :--- | :--- |")?;
    writeln!(
        w,
        "| ProkaDiff Engine Version | `{}` |",
        audit.provenance.prokadiff_version
    )?;
    writeln!(w, "| Git Commit SHA | `{}` |", audit.provenance.git_commit)?;
    if let Some(ref_sha) = &audit.provenance.reference_sha256 {
        writeln!(w, "| Reference SHA256 | `{}` |", ref_sha)?;
    }
    if let Some(bt2) = &audit.provenance.bowtie2_version {
        writeln!(w, "| Bowtie2 Version | `{}` |", bt2)?;
    }
    writeln!(
        w,
        "| Offtarget Engine Status | `{}` |",
        audit.provenance.offtarget_search_status
    )?;
    writeln!(
        w,
        "| CFD Scoring Gate | `{}` |",
        audit.provenance.cfd_scoring_status
    )?;
    writeln!(
        w,
        "| Bulge Search Status | `{}` |",
        audit.provenance.bulge_search_status
    )?;
    writeln!(
        w,
        "| Analysis Timestamp | `{}` |\n",
        audit.provenance.run_timestamp
    )?;
    Ok(())
}
