use crate::{EditorKind, RefContig};

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
    let spacer: Vec<u8> = spacer.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let pam: Vec<u8> = pam.bytes().map(|b| b.to_ascii_uppercase()).collect();
    if spacer.is_empty() || pam.is_empty() {
        return Vec::new();
    }
    let mut sites = Vec::new();
    for rec in refs {
        let seq = &rec.seq;
        match editor {
            EditorKind::Dsb => {}
            EditorKind::Cas9 => scan_cas9(
                rec.name.as_str(),
                seq,
                &spacer,
                &pam,
                max_mismatches,
                &mut sites,
            ),
            EditorKind::Cas12a => scan_cas12a(
                rec.name.as_str(),
                seq,
                &spacer,
                &pam,
                max_mismatches,
                &mut sites,
            ),
        }
    }
    sites
}

fn scan_cas9(
    seq_id: &str,
    seq: &[u8],
    spacer: &[u8],
    pam: &[u8],
    max_mm: u32,
    out: &mut Vec<HomologSite>,
) {
    let sp = spacer.len();
    let pn = pam.len();
    if seq.len() < sp + pn {
        return;
    }
    let pam_rc = revcomp_iupac(pam);
    let pam_s = String::from_utf8_lossy(pam).into_owned();
    // Plus: [protospacer][PAM]
    for p in 0..=seq.len() - pn {
        if !iupac_eq_slice(pam, &seq[p..p + pn]) {
            continue;
        }
        if p < sp {
            continue;
        }
        let mm = hamming(&seq[p - sp..p], spacer);
        if mm <= max_mm {
            out.push(HomologSite {
                seq_id: seq_id.to_string(),
                start: (p - sp) as u64 + 1,
                end: (p + pn) as u64,
                strand: '+',
                mismatches: mm,
                pam: pam_s.clone(),
            });
        }
    }
    // Minus: [revcomp PAM][revcomp protospacer] on the top strand.
    for p in 0..=seq.len() - pn {
        if !iupac_eq_slice(&pam_rc, &seq[p..p + pn]) {
            continue;
        }
        let sp0 = p + pn;
        if sp0 + sp > seq.len() {
            continue;
        }
        let proto = revcomp_dna(&seq[sp0..sp0 + sp]);
        let mm = hamming(&proto, spacer);
        if mm <= max_mm {
            out.push(HomologSite {
                seq_id: seq_id.to_string(),
                start: p as u64 + 1,
                end: (sp0 + sp) as u64,
                strand: '-',
                mismatches: mm,
                pam: pam_s.clone(),
            });
        }
    }
}

fn scan_cas12a(
    seq_id: &str,
    seq: &[u8],
    spacer: &[u8],
    pam: &[u8],
    max_mm: u32,
    out: &mut Vec<HomologSite>,
) {
    let sp = spacer.len();
    let pn = pam.len();
    if seq.len() < sp + pn {
        return;
    }
    let pam_rc = revcomp_iupac(pam);
    let pam_s = String::from_utf8_lossy(pam).into_owned();
    // Plus: [PAM][spacer]
    for p in 0..=seq.len() - pn {
        if !iupac_eq_slice(pam, &seq[p..p + pn]) {
            continue;
        }
        let sp0 = p + pn;
        if sp0 + sp > seq.len() {
            continue;
        }
        let mm = hamming(&seq[sp0..sp0 + sp], spacer);
        if mm <= max_mm {
            out.push(HomologSite {
                seq_id: seq_id.to_string(),
                start: p as u64 + 1,
                end: (sp0 + sp) as u64,
                strand: '+',
                mismatches: mm,
                pam: pam_s.clone(),
            });
        }
    }
    // Minus: [revcomp spacer][revcomp PAM] on the top strand.
    for p in 0..=seq.len() - pn {
        if !iupac_eq_slice(&pam_rc, &seq[p..p + pn]) {
            continue;
        }
        if p < sp {
            continue;
        }
        let proto = revcomp_dna(&seq[p - sp..p]);
        let mm = hamming(&proto, spacer);
        if mm <= max_mm {
            out.push(HomologSite {
                seq_id: seq_id.to_string(),
                start: (p - sp) as u64 + 1,
                end: (p + pn) as u64,
                strand: '-',
                mismatches: mm,
                pam: pam_s.clone(),
            });
        }
    }
}

fn hamming(a: &[u8], b: &[u8]) -> u32 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| u32::from(!x.eq_ignore_ascii_case(y)))
        .sum()
}

fn iupac_eq_slice(pat: &[u8], seq: &[u8]) -> bool {
    pat.len() == seq.len() && pat.iter().zip(seq.iter()).all(|(p, s)| iupac_eq(*p, *s))
}

fn iupac_eq(pat: u8, base: u8) -> bool {
    let b = base.to_ascii_uppercase();
    match pat.to_ascii_uppercase() {
        x if x == b && matches!(x, b'A' | b'C' | b'G' | b'T') => true,
        b'A' | b'C' | b'G' | b'T' => false,
        b'N' => matches!(b, b'A' | b'C' | b'G' | b'T'),
        b'V' => matches!(b, b'A' | b'C' | b'G'),
        b'B' => matches!(b, b'C' | b'G' | b'T'),
        b'D' => matches!(b, b'A' | b'G' | b'T'),
        b'H' => matches!(b, b'A' | b'C' | b'T'),
        b'R' => matches!(b, b'A' | b'G'),
        b'Y' => matches!(b, b'C' | b'T'),
        b'W' => matches!(b, b'A' | b'T'),
        b'S' => matches!(b, b'C' | b'G'),
        b'K' => matches!(b, b'G' | b'T'),
        b'M' => matches!(b, b'A' | b'C'),
        _ => false,
    }
}

fn complement_iupac(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        b'V' => b'B',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'R' => b'Y',
        b'Y' => b'R',
        b'W' => b'W',
        b'S' => b'S',
        b'K' => b'M',
        b'M' => b'K',
        x => x,
    }
}

fn revcomp_iupac(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().copied().map(complement_iupac).collect()
}

fn revcomp_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            x => x,
        })
        .collect()
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
