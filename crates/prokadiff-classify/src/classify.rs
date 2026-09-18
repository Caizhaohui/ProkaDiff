use prokadiff_gd::{GdEntry, GdKind, GenomeDiff};

use crate::homolog::{scan_homologs, HomologSite};
use crate::intended::{entry_intervals, mask_intended, IntendedEdit};
use crate::{is_product_mutation, is_structural, EditorKind, RefContig};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MutationClass {
    Structural,
    NearHomolog,
    ScatteredSnv,
}

impl MutationClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Structural => "structural",
            Self::NearHomolog => "near_homolog",
            Self::ScatteredSnv => "scattered_snv",
        }
    }
}

#[derive(Clone, Debug)]
pub struct ClassifyOptions {
    pub editor: EditorKind,
    pub spacer: Option<String>,
    pub pam: Option<String>,
    pub near_distance: u64,
    pub max_mismatches: u32,
    pub hypothesis: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClassifiedMutation {
    pub entry: GdEntry,
    pub class: MutationClass,
    pub pam_profile: Option<String>,
    pub offtarget_mismatch: Option<u32>,
    pub distance_to_site: Option<u64>,
    pub hypothesis: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct ClassifyResult {
    pub unintended: Vec<ClassifiedMutation>,
    pub intended_observed: Vec<GdEntry>,
    /// Number of rows in the `--intended` table (`intended.len()`). Zero when omitted.
    pub intended_declared: usize,
    pub starter_vs_ref: usize,
    /// Per-edit assessment (FIX-015): one entry per row in the `--intended` table.
    /// `None` when no intended table was provided.
    pub intended_edit_assessments: Option<Vec<crate::intended::IntendedEditAssessment>>,
}

pub fn classify(
    edited: &GenomeDiff,
    starter: &GenomeDiff,
    intended: &[IntendedEdit],
    refs: &[RefContig],
    opts: &ClassifyOptions,
) -> ClassifyResult {
    let starter_muts: Vec<GdEntry> = starter
        .entries
        .iter()
        .filter(|e| is_product_mutation(e.kind))
        .cloned()
        .collect();
    let edited_only = GenomeDiff {
        metadata: edited.metadata.clone(),
        entries: edited
            .entries
            .iter()
            .filter(|e| is_product_mutation(e.kind))
            .cloned()
            .collect(),
    };
    let starter_only = GenomeDiff {
        metadata: starter.metadata.clone(),
        entries: starter_muts.clone(),
    };
    let diff = edited_only.subtract(&starter_only);
    let (remain, observed) = mask_intended(&diff.entries, intended);

    let pam_used = resolved_pam(opts);
    let sites = match (opts.editor, opts.spacer.as_deref(), pam_used.as_deref()) {
        (EditorKind::Dsb, _, _) | (_, None | Some(""), _) | (_, _, None | Some("")) => Vec::new(),
        (editor, Some(spacer), Some(pam)) => {
            scan_homologs(refs, spacer, pam, editor, opts.max_mismatches)
        }
    };

    let mob_positions: Vec<(String, u64)> = remain
        .iter()
        .filter(|e| e.kind == GdKind::Mob)
        .filter_map(|e| Some((e.seq_id()?.to_string(), e.position()?)))
        .collect();

    let mut unintended = Vec::with_capacity(remain.len());
    for e in remain {
        if e.kind == GdKind::Jc && e.attrs.contains_key("mob_evidence") {
            if let (Some(seq), Some(pos)) = (e.seq_id(), e.position()) {
                if mob_positions
                    .iter()
                    .any(|(m_seq, m_pos)| m_seq == seq && pos.abs_diff(*m_pos) <= 20)
                {
                    continue;
                }
            }
        }
        unintended.push(label_one(e, &sites, opts, pam_used.as_deref()));
    }

    // FIX-015: build per-edit assessments using all mutations (including those in unintended)
    //
    // RW-004 fix: the intended-edit assessment must also see `MC` (missing-coverage)
    // evidence at declared loci, not just product mutations. `MC` is intentionally
    // excluded from `is_product_mutation` (it stays evidence, not a classified
    // mutation, so `unintended.tsv`/`summary.txt` are unaffected — see docs/schema.md),
    // but that meant a large aberrant on-target deletion with no matching on-target
    // JC/DEL call (e.g. an IS-mediated deletion where the on-target junction itself
    // does not survive JC accept thresholds) was invisible to `assess_intended_edits`,
    // and got reported as `Missing` instead of `UnexpectedStructure`. Only edited-strain
    // MC spans that overlap a declared intended-edit locus are added here — this does
    // not touch the general JC/MC genome-wide noise (tracked separately as RW-002/RW-003)
    // and does not change `diff`/`unintended` output for any input.
    let mut all_diff_entries: Vec<GdEntry> = diff.entries.clone();
    all_diff_entries.extend(
        edited
            .entries
            .iter()
            .filter(|e| e.kind == GdKind::Mc && mc_overlaps_any_intended(e, intended))
            .cloned(),
    );
    let intended_edit_assessments = if intended.is_empty() {
        None
    } else {
        use crate::intended::assess_intended_edits;
        Some(assess_intended_edits(&all_diff_entries, intended))
    };

    ClassifyResult {
        unintended,
        intended_observed: observed.into_iter().cloned().collect(),
        intended_declared: intended.len(),
        starter_vs_ref: starter_muts.len(),
        intended_edit_assessments,
    }
}

fn resolved_pam(opts: &ClassifyOptions) -> Option<String> {
    if let Some(p) = opts.pam.as_deref() {
        if !p.is_empty() {
            return Some(p.to_ascii_uppercase());
        }
    }
    opts.editor.default_pam().map(str::to_string)
}

fn label_one(
    e: &GdEntry,
    sites: &[HomologSite],
    opts: &ClassifyOptions,
    pam_used: Option<&str>,
) -> ClassifiedMutation {
    let del_size = if e.kind == GdKind::Del {
        e.fields.get(2).and_then(|s| s.parse().ok())
    } else {
        None
    };
    let (class, pam_profile, offtarget_mismatch, distance_to_site) =
        if is_structural(e.kind, del_size) {
            (MutationClass::Structural, None, None, None)
        } else if opts.editor != EditorKind::Dsb {
            if let Some((dist, mm)) = nearest_site(e, sites, opts.near_distance) {
                (
                    MutationClass::NearHomolog,
                    pam_used.map(str::to_string),
                    Some(mm),
                    Some(dist),
                )
            } else {
                (MutationClass::ScatteredSnv, None, None, None)
            }
        } else {
            (MutationClass::ScatteredSnv, None, None, None)
        };
    let hypothesis = if opts.hypothesis && class == MutationClass::ScatteredSnv {
        Some("sos_widney2014".into())
    } else {
        None
    };
    ClassifiedMutation {
        entry: e.clone(),
        class,
        pam_profile,
        offtarget_mismatch,
        distance_to_site,
        hypothesis,
    }
}

/// RW-004: whether an `MC` entry's interval overlaps any declared intended-edit locus.
/// Used only to widen the input to `assess_intended_edits`; does not affect `unintended`.
fn mc_overlaps_any_intended(e: &GdEntry, intended: &[IntendedEdit]) -> bool {
    entry_intervals(e).iter().any(|(sid, a, b)| {
        intended
            .iter()
            .any(|t| sid == &t.seq_id && *a <= t.end && t.start <= *b)
    })
}

fn nearest_site(e: &GdEntry, sites: &[HomologSite], max_d: u64) -> Option<(u64, u32)> {
    let mut best: Option<(u64, u32)> = None;
    for (sid, a, b) in entry_intervals(e) {
        for s in sites {
            if s.seq_id != sid {
                continue;
            }
            let d = interval_distance(a, b, s.start, s.end);
            if d > max_d {
                continue;
            }
            match best {
                None => best = Some((d, s.mismatches)),
                Some((bd, bm)) if d < bd || (d == bd && s.mismatches < bm) => {
                    best = Some((d, s.mismatches));
                }
                _ => {}
            }
        }
    }
    best
}

fn interval_distance(a0: u64, a1: u64, b0: u64, b1: u64) -> u64 {
    if a0 <= b1 && b0 <= a1 {
        0
    } else if a0 > b1 {
        a0 - b1
    } else {
        b0 - a1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intended::{parse_intended, IntendedEditStatus};
    use prokadiff_gd::GenomeDiff;

    fn gd(entries: Vec<GdEntry>) -> GenomeDiff {
        GenomeDiff {
            metadata: vec![("GENOME_DIFF".into(), "1.0".into())],
            entries,
        }
    }

    fn dsb_opts() -> ClassifyOptions {
        ClassifyOptions {
            editor: EditorKind::Dsb,
            spacer: None,
            pam: None,
            near_distance: 5000,
            max_mismatches: 3,
            hypothesis: false,
        }
    }

    /// RW-004 regression test: an on-target aberrant large deletion that shows up ONLY
    /// as an `MC` (missing-coverage) span in the edited sample -- no matching on-target
    /// JC/DEL call survives to the diff -- must still classify the declared `del` edit
    /// as `UnexpectedStructure`, not `Missing`, when run through the full `classify()`
    /// entrypoint (not just `assess_single_edit`/`assess_intended_edits` directly).
    ///
    /// This mirrors the real BL21 `B21_3_1` case: declared edit is the clean 860 bp
    /// `lacZ` deletion (334876-335735), but the observed event is a much larger
    /// ~20 kb missing-coverage span (331956-352202) with no on-target JC/DEL emitted.
    #[test]
    fn mc_only_aberrant_deletion_at_intended_locus_yields_unexpected_structure() {
        let edited = gd(vec![GdEntry::mc(903, "chr", 331956, 352202, 0, 0)]);
        let starter = gd(vec![]);
        let intended =
            parse_intended("seq_id\tstart\tend\tref\talt\tkind\nchr\t334876\t335735\t.\t.\tdel\n")
                .unwrap();

        let out = classify(&edited, &starter, &intended, &[], &dsb_opts());

        let assessments = out
            .intended_edit_assessments
            .expect("intended edits were declared");
        assert_eq!(assessments.len(), 1);
        assert_eq!(
            assessments[0].status,
            IntendedEditStatus::UnexpectedStructure
        );

        // The MC-widening must not leak into the general unintended/diff output: MC
        // stays evidence, not a classified mutation (unintended.tsv backward compat).
        assert!(out.unintended.is_empty());
    }

    /// An MC span far from any declared intended-edit locus must NOT be pulled into
    /// the intended-edit assessment (the widening in `classify()` is locus-scoped).
    #[test]
    fn mc_far_from_intended_locus_does_not_affect_assessment() {
        let edited = gd(vec![GdEntry::mc(1, "chr", 900_000, 950_000, 0, 0)]);
        let starter = gd(vec![]);
        let intended =
            parse_intended("seq_id\tstart\tend\tref\talt\tkind\nchr\t334876\t335735\t.\t.\tdel\n")
                .unwrap();

        let out = classify(&edited, &starter, &intended, &[], &dsb_opts());

        let assessments = out
            .intended_edit_assessments
            .expect("intended edits were declared");
        assert_eq!(assessments.len(), 1);
        assert_eq!(assessments[0].status, IntendedEditStatus::Missing);
    }
}
