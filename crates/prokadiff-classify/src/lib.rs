//! Product-layer: `M_edited \ M_parent \ {intended}`, then classes (3)→(1)→(2).

#![deny(unsafe_code)]

mod classify;
mod homolog;
mod intended;

pub use classify::{classify, ClassifiedMutation, ClassifyOptions, ClassifyResult, MutationClass};
pub use homolog::{scan_homologs, HomologSite, DEFAULT_MAX_MISMATCHES, DEFAULT_NEAR_DISTANCE};
pub use intended::{
    mask_intended, parse_intended, parse_intended_path, IntendedEdit, IntendedError,
};

use prokadiff_gd::GdKind;

/// First-period editor. CAST / IS110 are rejected in the CLI, not here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EditorKind {
    Cas9,
    Cas12a,
    Dsb,
}

impl EditorKind {
    pub fn default_pam(self) -> Option<&'static str> {
        match self {
            Self::Cas9 => Some("NGG"),
            Self::Cas12a => Some("TTTV"),
            Self::Dsb => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Cas9 => "cas9",
            Self::Cas12a => "cas12a",
            Self::Dsb => "dsb",
        }
    }
}

/// One reference contig (coordinate skeleton). Classify does not depend on BAM.
#[derive(Clone, Debug)]
pub struct RefContig {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Mutation types that enter the unintended list (not RA/MC/UN evidence).
pub fn is_product_mutation(kind: GdKind) -> bool {
    matches!(
        kind,
        GdKind::Snp
            | GdKind::Ins
            | GdKind::Del
            | GdKind::Mob
            | GdKind::Amp
            | GdKind::Con
            | GdKind::Jc
    )
}

/// Structural (class 3): MOB/JC/AMP/CON, or DEL longer than the RA 2 bp subset.
pub fn is_structural(kind: GdKind, del_size: Option<u64>) -> bool {
    match kind {
        GdKind::Mob | GdKind::Jc | GdKind::Amp | GdKind::Con => true,
        GdKind::Del => del_size.is_some_and(|s| s > 2),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use prokadiff_gd::{GdEntry, GenomeDiff};

    fn gd(entries: Vec<GdEntry>) -> GenomeDiff {
        GenomeDiff {
            metadata: vec![("GENOME_DIFF".into(), "1.0".into())],
            entries,
        }
    }

    fn cas9_opts(hypothesis: bool) -> ClassifyOptions {
        ClassifyOptions {
            editor: EditorKind::Cas9,
            spacer: Some("GACTGACTGACTGACTGACT".into()),
            pam: None,
            near_distance: DEFAULT_NEAR_DISTANCE,
            max_mismatches: DEFAULT_MAX_MISMATCHES,
            hypothesis,
        }
    }

    /// 10 bp prefix + 20 bp spacer + TGG + poly-A.
    fn cas9_ref() -> Vec<RefContig> {
        let spacer = b"GACTGACTGACTGACTGACT";
        let mut seq = b"AAAAAAAAAA".to_vec();
        seq.extend_from_slice(spacer);
        seq.extend_from_slice(b"TGG");
        seq.extend(std::iter::repeat_n(b'A', 200));
        vec![RefContig {
            name: "chr".into(),
            seq,
        }]
    }

    #[test]
    fn fifty_historical_snps_are_subtracted() {
        let hist: Vec<GdEntry> = (1..=50)
            .map(|i| GdEntry::snp(i, "chr", u64::from(i) * 10, "A"))
            .collect();
        let mut edited_e = hist.clone();
        edited_e.push(GdEntry::snp(99, "chr", 8000, "T"));
        let out = classify(&gd(edited_e), &gd(hist), &[], &[], &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1);
        assert_eq!(out.unintended[0].entry.fields[1], "8000");
        assert_eq!(out.unintended[0].class, MutationClass::ScatteredSnv);
        assert_eq!(out.starter_vs_ref, 50);
    }

    #[test]
    fn intended_snp_is_masked_not_unintended() {
        let edited = gd(vec![
            GdEntry::snp(1, "chr", 100, "G"),
            GdEntry::snp(2, "chr", 200, "C"),
        ]);
        let intended =
            parse_intended("seq_id\tstart\tend\tref\talt\tkind\nchr\t100\t100\tA\tG\tsnp\n")
                .unwrap();
        let out = classify(&edited, &gd(vec![]), &intended, &[], &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1);
        assert_eq!(out.unintended[0].entry.fields[1], "200");
        assert_eq!(out.intended_observed.len(), 1);
        assert_eq!(out.intended_observed[0].fields[1], "100");
        assert_eq!(out.intended_declared, 1);
    }

    #[test]
    fn classify_sets_intended_declared_to_slice_len_even_if_partial() {
        let edited = gd(vec![GdEntry::snp(1, "chr", 100, "G")]);
        let intended = parse_intended(
            "seq_id\tstart\tend\tref\talt\tkind\nchr\t100\t100\tA\tG\tsnp\nchr\t200\t200\tA\tC\tsnp\n",
        )
        .unwrap();
        let out = classify(&edited, &gd(vec![]), &intended, &[], &cas9_opts(false));
        assert_eq!(out.intended_declared, 2);
        assert_eq!(out.intended_observed.len(), 1);
        assert_eq!(out.unintended.len(), 0);
    }

    #[test]
    fn omitted_intended_keeps_all_post_subtract() {
        let edited = gd(vec![GdEntry::snp(1, "chr", 100, "G")]);
        let out = classify(&edited, &gd(vec![]), &[], &[], &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1);
        assert!(out.intended_observed.is_empty());
        assert_eq!(out.intended_declared, 0);
    }

    #[test]
    fn mob_is_structural_even_when_near_cas9_site() {
        let refs = cas9_ref();
        let edited = gd(vec![GdEntry::mob(1, "chr", 20, "IS1", "+", 8)]);
        let out = classify(&edited, &gd(vec![]), &[], &refs, &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1);
        assert_eq!(out.unintended[0].class, MutationClass::Structural);
    }

    #[test]
    fn large_del_is_structural_short_indel_is_not() {
        let edited = gd(vec![
            GdEntry::del(1, "chr", 10, 2),
            GdEntry::del(2, "chr", 80, 40),
            GdEntry::ins(3, "chr", 5, "AT"),
        ]);
        let out = classify(&edited, &gd(vec![]), &[], &[], &cas9_opts(false));
        let classes: Vec<_> = out
            .unintended
            .iter()
            .map(|c| (c.entry.kind, c.class))
            .collect();
        assert!(classes.contains(&(GdKind::Del, MutationClass::ScatteredSnv)));
        assert!(classes.contains(&(GdKind::Del, MutationClass::Structural)));
        assert!(classes.contains(&(GdKind::Ins, MutationClass::ScatteredSnv)));
    }

    #[test]
    fn cas9_snp_within_50bp_of_ngg_homolog_is_near_homolog() {
        let refs = cas9_ref();
        // Site is 1-based 11..=33; position 40 is 7 bp past PAM.
        let edited = gd(vec![GdEntry::snp(1, "chr", 40, "C")]);
        let out = classify(&edited, &gd(vec![]), &[], &refs, &cas9_opts(false));
        assert_eq!(out.unintended[0].class, MutationClass::NearHomolog);
        assert_eq!(out.unintended[0].pam_profile.as_deref(), Some("NGG"));
        assert_eq!(out.unintended[0].offtarget_mismatch, Some(0));
        assert_eq!(out.unintended[0].distance_to_site, Some(7));
    }

    #[test]
    fn dsb_does_not_label_near_homolog_but_still_reports_scattered() {
        let refs = cas9_ref();
        let edited = gd(vec![GdEntry::snp(1, "chr", 40, "C")]);
        let mut opts = cas9_opts(false);
        opts.editor = EditorKind::Dsb;
        opts.spacer = None;
        let out = classify(&edited, &gd(vec![]), &[], &refs, &opts);
        assert_eq!(out.unintended[0].class, MutationClass::ScatteredSnv);
        assert!(out.unintended[0].pam_profile.is_none());
    }

    #[test]
    fn cas9_snp_far_from_site_is_scattered() {
        let refs = cas9_ref();
        let edited = gd(vec![GdEntry::snp(1, "chr", 200, "C")]);
        let out = classify(&edited, &gd(vec![]), &[], &refs, &cas9_opts(false));
        assert_eq!(out.unintended[0].class, MutationClass::ScatteredSnv);
        assert!(out.unintended[0].distance_to_site.is_none());
    }

    #[test]
    fn cas12a_tttv_homolog_labels_near() {
        let spacer = b"GACTGACTGACTGACTGACT";
        let mut seq = b"AAAAAAAAAA".to_vec();
        seq.extend_from_slice(b"TTTA");
        seq.extend_from_slice(spacer);
        seq.extend(std::iter::repeat_n(b'A', 80));
        let refs = vec![RefContig {
            name: "chr".into(),
            seq,
        }];
        let edited = gd(vec![GdEntry::snp(1, "chr", 40, "C")]);
        let opts = ClassifyOptions {
            editor: EditorKind::Cas12a,
            spacer: Some("GACTGACTGACTGACTGACT".into()),
            pam: None,
            near_distance: DEFAULT_NEAR_DISTANCE,
            max_mismatches: DEFAULT_MAX_MISMATCHES,
            hypothesis: false,
        };
        let out = classify(&edited, &gd(vec![]), &[], &refs, &opts);
        assert_eq!(out.unintended[0].class, MutationClass::NearHomolog);
        assert_eq!(out.unintended[0].pam_profile.as_deref(), Some("TTTV"));
    }

    #[test]
    fn scattered_gets_sos_hypothesis_unless_disabled() {
        let edited = gd(vec![GdEntry::snp(1, "chr", 200, "C")]);
        let with_h = classify(&edited, &gd(vec![]), &[], &cas9_ref(), &cas9_opts(true));
        assert_eq!(
            with_h.unintended[0].hypothesis.as_deref(),
            Some("sos_widney2014")
        );
        let no_h = classify(&edited, &gd(vec![]), &[], &cas9_ref(), &cas9_opts(false));
        assert!(no_h.unintended[0].hypothesis.is_none());
    }

    #[test]
    fn ra_mc_un_are_not_unintended_mutations() {
        let edited = gd(vec![
            GdEntry::snp(1, "chr", 50, "A"),
            GdEntry::mc(2, "chr", 1, 10, 0, 0),
            GdEntry {
                kind: GdKind::Un,
                id: 3,
                parent_ids: vec![],
                fields: vec!["chr".into(), "1".into(), "5".into()],
                attrs: Default::default(),
            },
            GdEntry {
                kind: GdKind::Ra,
                id: 4,
                parent_ids: vec![],
                fields: vec![
                    "chr".into(),
                    "50".into(),
                    "0".into(),
                    "T".into(),
                    "A".into(),
                ],
                attrs: Default::default(),
            },
        ]);
        let out = classify(&edited, &gd(vec![]), &[], &[], &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1);
        assert_eq!(out.unintended[0].entry.kind, GdKind::Snp);
    }

    #[test]
    fn cassette_intended_masks_overlapping_jc() {
        let edited = gd(vec![GdEntry::jc(1, "chr", 100, "+", "chr", 500, "-", 0)]);
        let intended =
            parse_intended("seq_id\tstart\tend\tref\talt\tkind\nchr\t90\t110\t.\t.\tcassette\n")
                .unwrap();
        let out = classify(&edited, &gd(vec![]), &intended, &[], &cas9_opts(false));
        assert!(out.unintended.is_empty());
        assert_eq!(out.intended_observed.len(), 1);
    }

    #[test]
    fn intended_del_does_not_mask_insertion() {
        // User declared a 2 bp DEL; an INS in the same interval must stay unintended.
        let edited = gd(vec![GdEntry::ins(1, "chr", 100, "AT")]);
        let intended =
            parse_intended("seq_id\tstart\tend\tref\talt\tkind\nchr\t100\t101\t.\t.\tdel\n")
                .unwrap();
        let out = classify(&edited, &gd(vec![]), &intended, &[], &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1);
        assert!(out.intended_observed.is_empty());
    }

    #[test]
    fn intended_del_requires_matching_span() {
        let edited = gd(vec![GdEntry::del(1, "chr", 100, 5)]);
        let declared = |span_end: u64| {
            parse_intended(&format!(
                "seq_id\tstart\tend\tref\talt\tkind\nchr\t100\t{span_end}\t.\t.\tdel\n"
            ))
            .unwrap()
        };
        let out = classify(&edited, &gd(vec![]), &declared(102), &[], &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1, "2 bp declared vs 5 bp observed");
        assert!(out.intended_observed.is_empty());
        let out = classify(&edited, &gd(vec![]), &declared(104), &[], &cas9_opts(false));
        assert!(out.unintended.is_empty());
        assert_eq!(out.intended_observed.len(), 1);
    }

    #[test]
    fn intended_del_with_wrong_alt_does_not_mask() {
        let edited = gd(vec![GdEntry::del(1, "chr", 100, 2)]);
        let intended =
            parse_intended("seq_id\tstart\tend\tref\talt\tkind\nchr\t100\t101\t.\tG\tdel\n")
                .unwrap();
        let out = classify(&edited, &gd(vec![]), &intended, &[], &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1);
        assert!(out.intended_observed.is_empty());
    }

    #[test]
    fn unknown_intended_kind_does_not_mask_snp() {
        let edited = gd(vec![GdEntry::snp(1, "chr", 100, "G")]);
        let intended =
            parse_intended("seq_id\tstart\tend\tref\talt\tkind\nchr\t100\t100\tA\tG\tweird\n")
                .unwrap();
        let out = classify(&edited, &gd(vec![]), &intended, &[], &cas9_opts(false));
        assert_eq!(out.unintended.len(), 1);
        assert!(out.intended_observed.is_empty());
    }

    #[test]
    fn cas9_minus_strand_snp_is_near_homolog() {
        let spacer_rc = b"AGTCAGTCAGTCAGTCAGTC";
        let mut seq = b"AAAAAAAAAA".to_vec();
        seq.extend_from_slice(b"CCA");
        seq.extend_from_slice(spacer_rc);
        seq.extend(std::iter::repeat_n(b'A', 200));
        let refs = vec![RefContig {
            name: "chr".into(),
            seq,
        }];
        let edited = gd(vec![GdEntry::snp(1, "chr", 40, "C")]);
        let out = classify(&edited, &gd(vec![]), &[], &refs, &cas9_opts(false));
        assert_eq!(out.unintended[0].class, MutationClass::NearHomolog);
        assert_eq!(out.unintended[0].pam_profile.as_deref(), Some("NGG"));
        assert_eq!(out.unintended[0].offtarget_mismatch, Some(0));
        assert_eq!(out.unintended[0].distance_to_site, Some(7));
    }

    #[test]
    fn distance_50bp_is_near_homolog_51bp_is_scattered() {
        // cas9_ref site occupies 11..=33; 33+50=83, 33+51=84.
        let edited_50 = gd(vec![GdEntry::snp(1, "chr", 83, "C")]);
        let out = classify(&edited_50, &gd(vec![]), &[], &cas9_ref(), &cas9_opts(false));
        assert_eq!(out.unintended[0].class, MutationClass::NearHomolog);
        assert_eq!(out.unintended[0].distance_to_site, Some(50));
        let edited_51 = gd(vec![GdEntry::snp(1, "chr", 84, "C")]);
        let out = classify(&edited_51, &gd(vec![]), &[], &cas9_ref(), &cas9_opts(false));
        assert_eq!(out.unintended[0].class, MutationClass::ScatteredSnv);
    }

    #[test]
    fn pam_override_nag_matches_tag_not_default_ngg() {
        let spacer = b"GACTGACTGACTGACTGACT";
        let mut seq = b"AAAAAAAAAA".to_vec();
        seq.extend_from_slice(spacer);
        seq.extend_from_slice(b"TAG");
        seq.extend(std::iter::repeat_n(b'A', 200));
        let refs = vec![RefContig {
            name: "chr".into(),
            seq,
        }];
        let edited = gd(vec![GdEntry::snp(1, "chr", 40, "C")]);
        let default_ngg = classify(&edited, &gd(vec![]), &[], &refs, &cas9_opts(false));
        assert_eq!(default_ngg.unintended[0].class, MutationClass::ScatteredSnv);
        let mut nag = cas9_opts(false);
        nag.pam = Some("NAG".into());
        let out = classify(&edited, &gd(vec![]), &[], &refs, &nag);
        assert_eq!(out.unintended[0].class, MutationClass::NearHomolog);
        assert_eq!(out.unintended[0].pam_profile.as_deref(), Some("NAG"));
    }
}
