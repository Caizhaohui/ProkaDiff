use crate::{EditorKind, RefContig};
use prokadiff_offtarget::model::{NucleaseProfile, PamSide};
use prokadiff_offtarget::rust_search::scan_contig;

pub const DEFAULT_NEAR_DISTANCE: u64 = 50;
pub const DEFAULT_MAX_MISMATCHES: u32 = 4;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HomologSite {
    pub seq_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub mismatches: u32,
    pub pam: String,
}

pub fn scan_homologs(
    refs: &[RefContig],
    spacer: &str,
    pam: &str,
    editor: EditorKind,
    max_mismatches: u32,
) -> Vec<HomologSite> {
    if spacer.trim().is_empty() || pam.trim().is_empty() {
        return Vec::new();
    }
    let profile = match editor {
        EditorKind::Dsb => return Vec::new(),
        EditorKind::Cas9 => NucleaseProfile::custom("Cas9", spacer.len(), pam, PamSide::ThreePrime),
        EditorKind::Cas12a => {
            NucleaseProfile::custom("Cas12a", spacer.len(), pam, PamSide::FivePrime)
        }
    };

    let mut sites = Vec::new();
    for rec in refs {
        let found = scan_contig(
            &rec.name,
            &rec.seq,
            spacer,
            &profile,
            max_mismatches,
            sites.len(),
        );
        for s in found {
            sites.push(HomologSite {
                seq_id: s.seq_id,
                start: s.start,
                end: s.end,
                strand: s.strand.as_char(),
                mismatches: s.mismatches,
                pam: s.pam,
            });
        }
    }
    sites
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RefContig;

    const SPACER: &str = "GACTGACTGACTGACTGACT";
    /// Reverse complement of `SPACER`.
    const SPACER_RC: &[u8] = b"AGTCAGTCAGTCAGTCAGTC";

    fn contig(seq: Vec<u8>) -> Vec<RefContig> {
        vec![RefContig {
            name: "chr".into(),
            seq,
        }]
    }

    #[test]
    fn cas9_minus_strand_ngg_is_scanned() {
        // Top strand: [revcomp PAM=CCA][revcomp spacer]  ≡ minus-strand spacer+TGG.
        let mut seq = b"AAAAAAAAAA".to_vec();
        seq.extend_from_slice(b"CCA");
        seq.extend_from_slice(SPACER_RC);
        seq.extend(std::iter::repeat_n(b'A', 80));
        let sites = scan_homologs(&contig(seq), SPACER, "NGG", EditorKind::Cas9, 0);
        assert_eq!(sites.len(), 1, "{sites:?}");
        assert_eq!(sites[0].strand, '-');
        assert_eq!(sites[0].mismatches, 0);
        assert_eq!(sites[0].start, 11);
        assert_eq!(sites[0].end, 33);
    }

    #[test]
    fn cas12a_minus_strand_tttv_is_scanned() {
        // Top strand: [revcomp spacer][revcomp TTTA=TAAA].
        let mut seq = b"AAAAAAAAAA".to_vec();
        seq.extend_from_slice(SPACER_RC);
        seq.extend_from_slice(b"TAAA");
        seq.extend(std::iter::repeat_n(b'A', 80));
        let sites = scan_homologs(&contig(seq), SPACER, "TTTV", EditorKind::Cas12a, 0);
        assert_eq!(sites.len(), 1, "{sites:?}");
        assert_eq!(sites[0].strand, '-');
        assert_eq!(sites[0].mismatches, 0);
        assert_eq!(sites[0].start, 11);
        assert_eq!(sites[0].end, 34);
    }

    #[test]
    fn four_spacer_mismatches_accepted_five_rejected() {
        let mut seq = b"AAAAAAAAAA".to_vec();
        seq.extend_from_slice(b"TCAGGACTGACTGACTGACT"); // GACT→TCAG: 4 real mismatches
        seq.extend_from_slice(b"TGG");
        seq.extend(std::iter::repeat_n(b'A', 80));
        let four = scan_homologs(&contig(seq.clone()), SPACER, "NGG", EditorKind::Cas9, 4);
        assert_eq!(four.len(), 1);
        assert_eq!(four[0].mismatches, 4);
        assert!(scan_homologs(&contig(seq.clone()), SPACER, "NGG", EditorKind::Cas9, 3).is_empty());

        seq[10..15].copy_from_slice(b"TCAGT"); // 5 mismatches (4 + spacer[4] G→T)
        let five = scan_homologs(&contig(seq), SPACER, "NGG", EditorKind::Cas9, 4);
        assert!(five.is_empty(), "{five:?}");
    }
}
