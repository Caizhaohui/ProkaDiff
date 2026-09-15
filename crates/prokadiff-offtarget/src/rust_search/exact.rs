use crate::model::{BulgeType, NucleaseProfile, OffTargetSite, PamSide, Strand};
use crate::normalize::{
    clean_nucleotide_string, hamming, iupac_match_slice, revcomp_dna, revcomp_iupac,
};

/// Scans reference contigs for candidate off-target sites using exact mismatch counting.
///
/// Fully replicates Cas-OFFinder mismatch-only search behavior across both strands and IUPAC PAMs.
pub fn scan_genome(
    refs: &[(String, Vec<u8>)],
    guide: &str,
    profile: &NucleaseProfile,
    max_mismatches: u32,
) -> Vec<OffTargetSite> {
    let mut all_sites = Vec::new();
    for (name, seq) in refs {
        let sites = scan_contig(name, seq, guide, profile, max_mismatches, all_sites.len());
        all_sites.extend(sites);
    }
    all_sites
}

/// Scans a single contig for candidate sites matching the nuclease profile and guide.
pub fn scan_contig(
    seq_id: &str,
    seq: &[u8],
    guide: &str,
    profile: &NucleaseProfile,
    max_mismatches: u32,
    id_offset: usize,
) -> Vec<OffTargetSite> {
    let guide_clean = clean_nucleotide_string(guide);
    let guide_bytes = guide_clean.as_bytes();
    let pam_pattern = clean_nucleotide_string(&profile.pam_pattern);
    let pam_bytes = pam_pattern.as_bytes();

    let sp = guide_bytes.len();
    let pn = pam_bytes.len();
    if sp == 0 || pn == 0 || seq.len() < sp + pn {
        return Vec::new();
    }

    let pam_rc = revcomp_iupac(pam_bytes);
    let mut sites = Vec::new();

    match profile.pam_side {
        PamSide::ThreePrime => {
            // SpCas9-style: 5'-[protospacer]-[PAM]-3'
            // Plus strand: top strand is [protospacer][PAM]
            for p in sp..=seq.len() - pn {
                let pam_slice = &seq[p..p + pn];
                if !iupac_match_slice(pam_bytes, pam_slice) {
                    continue;
                }
                let proto_slice = &seq[p - sp..p];
                let mm = hamming(proto_slice, guide_bytes);
                if mm <= max_mismatches {
                    let site_idx = id_offset + sites.len() + 1;
                    let proto_str = String::from_utf8_lossy(proto_slice);
                    let pam_str = String::from_utf8_lossy(pam_slice);
                    let (cfd_score, hsu_score) = if sp == 20
                        && (profile.name.eq_ignore_ascii_case("spcas9")
                            || profile.name.eq_ignore_ascii_case("cas9"))
                    {
                        (
                            crate::scoring::calculate_cfd_score(&guide_clean, &proto_str, &pam_str),
                            crate::scoring::calculate_hsu_score(&guide_clean, &proto_str),
                        )
                    } else {
                        (None, None)
                    };
                    sites.push(OffTargetSite {
                        site_id: format!("site_{site_idx}"),
                        seq_id: seq_id.to_string(),
                        start: (p - sp) as u64 + 1,
                        end: (p + pn) as u64,
                        strand: Strand::Plus,
                        guide: guide_clean.clone(),
                        target_seq: String::from_utf8_lossy(&seq[p - sp..p + pn]).into_owned(),
                        pam: pam_str.into_owned(),
                        mismatches: mm,
                        bulge_type: BulgeType::None,
                        bulge_size: 0,
                        search_backend: "rust_exact".to_string(),
                        cfd_score,
                        hsu_score,
                    });
                }
            }

            // Minus strand: top strand is [revcomp(PAM)][revcomp(protospacer)]
            for p in 0..=seq.len().saturating_sub(pn + sp) {
                let pam_slice = &seq[p..p + pn];
                if !iupac_match_slice(&pam_rc, pam_slice) {
                    continue;
                }
                let proto_rc_slice = &seq[p + pn..p + pn + sp];
                let proto = revcomp_dna(proto_rc_slice);
                let mm = hamming(&proto, guide_bytes);
                if mm <= max_mismatches {
                    let site_idx = id_offset + sites.len() + 1;
                    let proto_str = String::from_utf8_lossy(&proto);
                    let pam_rc_bytes = revcomp_dna(pam_slice);
                    let pam_str = String::from_utf8_lossy(&pam_rc_bytes);
                    let (cfd_score, hsu_score) = if sp == 20
                        && (profile.name.eq_ignore_ascii_case("spcas9")
                            || profile.name.eq_ignore_ascii_case("cas9"))
                    {
                        (
                            crate::scoring::calculate_cfd_score(&guide_clean, &proto_str, &pam_str),
                            crate::scoring::calculate_hsu_score(&guide_clean, &proto_str),
                        )
                    } else {
                        (None, None)
                    };
                    sites.push(OffTargetSite {
                        site_id: format!("site_{site_idx}"),
                        seq_id: seq_id.to_string(),
                        start: p as u64 + 1,
                        end: (p + pn + sp) as u64,
                        strand: Strand::Minus,
                        guide: guide_clean.clone(),
                        target_seq: String::from_utf8_lossy(&seq[p..p + pn + sp]).into_owned(),
                        pam: pam_str.into_owned(),
                        mismatches: mm,
                        bulge_type: BulgeType::None,
                        bulge_size: 0,
                        search_backend: "rust_exact".to_string(),
                        cfd_score,
                        hsu_score,
                    });
                }
            }
        }
        PamSide::FivePrime => {
            // Cas12a-style: 5'-[PAM]-[protospacer]-3'
            // Plus strand: top strand is [PAM][protospacer]
            for p in 0..=seq.len().saturating_sub(pn + sp) {
                let pam_slice = &seq[p..p + pn];
                if !iupac_match_slice(pam_bytes, pam_slice) {
                    continue;
                }
                let proto_slice = &seq[p + pn..p + pn + sp];
                let mm = hamming(proto_slice, guide_bytes);
                if mm <= max_mismatches {
                    let site_idx = id_offset + sites.len() + 1;
                    sites.push(OffTargetSite {
                        site_id: format!("site_{site_idx}"),
                        seq_id: seq_id.to_string(),
                        start: p as u64 + 1,
                        end: (p + pn + sp) as u64,
                        strand: Strand::Plus,
                        guide: guide_clean.clone(),
                        target_seq: String::from_utf8_lossy(&seq[p..p + pn + sp]).into_owned(),
                        pam: String::from_utf8_lossy(pam_slice).into_owned(),
                        mismatches: mm,
                        bulge_type: BulgeType::None,
                        bulge_size: 0,
                        search_backend: "rust_exact".to_string(),
                        cfd_score: None,
                        hsu_score: None,
                    });
                }
            }

            // Minus strand: top strand is [revcomp(protospacer)][revcomp(PAM)]
            for p in sp..=seq.len() - pn {
                let pam_slice = &seq[p..p + pn];
                if !iupac_match_slice(&pam_rc, pam_slice) {
                    continue;
                }
                let proto_rc_slice = &seq[p - sp..p];
                let proto = revcomp_dna(proto_rc_slice);
                let mm = hamming(&proto, guide_bytes);
                if mm <= max_mismatches {
                    let site_idx = id_offset + sites.len() + 1;
                    sites.push(OffTargetSite {
                        site_id: format!("site_{site_idx}"),
                        seq_id: seq_id.to_string(),
                        start: (p - sp) as u64 + 1,
                        end: (p + pn) as u64,
                        strand: Strand::Minus,
                        guide: guide_clean.clone(),
                        target_seq: String::from_utf8_lossy(&seq[p - sp..p + pn]).into_owned(),
                        pam: String::from_utf8_lossy(&revcomp_dna(pam_slice)).into_owned(),
                        mismatches: mm,
                        bulge_type: BulgeType::None,
                        bulge_size: 0,
                        search_backend: "rust_exact".to_string(),
                        cfd_score: None,
                        hsu_score: None,
                    });
                }
            }
        }
    }

    sites
}
