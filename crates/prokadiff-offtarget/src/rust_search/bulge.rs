use crate::model::{BulgeType, NucleaseProfile, OffTargetSite, PamSide, Strand};
use crate::normalize::{clean_nucleotide_string, iupac_match_slice, revcomp_dna, revcomp_iupac};

/// Options controlling off-target search with optional DNA/RNA bulges.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BulgeSearchOptions {
    pub max_mismatches: u32,
    pub max_dna_bulge: u32,
    pub max_rna_bulge: u32,
}

impl Default for BulgeSearchOptions {
    fn default() -> Self {
        Self {
            max_mismatches: 3,
            max_dna_bulge: 0,
            max_rna_bulge: 0,
        }
    }
}

/// Computes the minimum number of mismatches between target and guide with a 1-bp DNA bulge
/// (one nucleotide inserted in target DNA relative to guide RNA).
/// `target` has length L + 1, `guide` has length L.
pub fn min_mismatches_dna_bulge_1(target: &[u8], guide: &[u8]) -> u32 {
    let l = guide.len();
    if target.len() != l + 1 || l == 0 {
        return u32::MAX;
    }

    let mut pref = vec![0u32; l + 1];
    for i in 0..l {
        pref[i + 1] = pref[i]
            + if target[i].eq_ignore_ascii_case(&guide[i]) {
                0
            } else {
                1
            };
    }

    let mut suff = vec![0u32; l + 1];
    for i in (0..l).rev() {
        suff[i] = suff[i + 1]
            + if target[i + 1].eq_ignore_ascii_case(&guide[i]) {
                0
            } else {
                1
            };
    }

    let mut min_mm = u32::MAX;
    for k in 0..=l {
        let mm = pref[k] + suff[k];
        if mm < min_mm {
            min_mm = mm;
        }
    }
    min_mm
}

/// Computes the minimum number of mismatches between target and guide with a 2-bp DNA bulge
/// (two nucleotides inserted in target DNA relative to guide RNA).
/// `target` has length L + 2, `guide` has length L.
pub fn min_mismatches_dna_bulge_2(target: &[u8], guide: &[u8]) -> u32 {
    let l = guide.len();
    if target.len() != l + 2 || l == 0 {
        return u32::MAX;
    }

    let mut min_mm = u32::MAX;
    for k1 in 0..=l {
        for k2 in k1 + 1..=l + 1 {
            let mut mm = 0u32;
            let mut g_idx = 0;
            for (t_idx, &t_byte) in target.iter().enumerate() {
                if t_idx == k1 || t_idx == k2 {
                    continue;
                }
                if g_idx < l && !t_byte.eq_ignore_ascii_case(&guide[g_idx]) {
                    mm += 1;
                }
                g_idx += 1;
            }
            if mm < min_mm {
                min_mm = mm;
            }
        }
    }
    min_mm
}

/// Computes the minimum number of mismatches between target and guide with a 1-bp RNA bulge
/// (one nucleotide deleted in target DNA relative to guide RNA / inserted in guide).
/// `target` has length L - 1, `guide` has length L.
pub fn min_mismatches_rna_bulge_1(target: &[u8], guide: &[u8]) -> u32 {
    let l = guide.len();
    if target.len() + 1 != l || l <= 1 {
        return u32::MAX;
    }

    let mut pref = vec![0u32; l];
    for i in 0..l - 1 {
        pref[i + 1] = pref[i]
            + if target[i].eq_ignore_ascii_case(&guide[i]) {
                0
            } else {
                1
            };
    }

    let mut suff = vec![0u32; l];
    for i in (0..l - 1).rev() {
        suff[i] = suff[i + 1]
            + if target[i].eq_ignore_ascii_case(&guide[i + 1]) {
                0
            } else {
                1
            };
    }

    let mut min_mm = u32::MAX;
    for k in 0..l {
        let mm = pref[k] + suff[k];
        if mm < min_mm {
            min_mm = mm;
        }
    }
    min_mm
}

/// Computes the minimum number of mismatches between target and guide with a 2-bp RNA bulge
/// (two nucleotides deleted in target DNA relative to guide RNA / inserted in guide).
/// `target` has length L - 2, `guide` has length L.
pub fn min_mismatches_rna_bulge_2(target: &[u8], guide: &[u8]) -> u32 {
    let l = guide.len();
    if target.len() + 2 != l || l <= 2 {
        return u32::MAX;
    }

    let mut min_mm = u32::MAX;
    for k1 in 0..l {
        for k2 in k1 + 1..l {
            let mut mm = 0u32;
            let mut t_idx = 0;
            for (g_idx, &g_byte) in guide.iter().enumerate() {
                if g_idx == k1 || g_idx == k2 {
                    continue;
                }
                if t_idx < target.len() && !target[t_idx].eq_ignore_ascii_case(&g_byte) {
                    mm += 1;
                }
                t_idx += 1;
            }
            if mm < min_mm {
                min_mm = mm;
            }
        }
    }
    min_mm
}

/// Scans reference contigs for off-target sites including mismatch and DNA/RNA bulge search.
pub fn scan_genome_bulge(
    refs: &[(String, Vec<u8>)],
    guide: &str,
    profile: &NucleaseProfile,
    opts: BulgeSearchOptions,
) -> Vec<OffTargetSite> {
    let mut all_sites = Vec::new();
    for (name, seq) in refs {
        let sites = scan_contig_bulge(name, seq, guide, profile, opts, all_sites.len());
        all_sites.extend(sites);
    }
    all_sites
}

/// Scans a single contig for candidate sites matching the nuclease profile and guide,
/// considering exact mismatch-only as well as DNA/RNA bulges up to specified limits.
pub fn scan_contig_bulge(
    seq_id: &str,
    seq: &[u8],
    guide: &str,
    profile: &NucleaseProfile,
    opts: BulgeSearchOptions,
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

    // 1. Run 0-bulge (exact mismatch) scan first.
    let exact_sites = crate::rust_search::exact::scan_contig(
        seq_id,
        seq,
        &guide_clean,
        profile,
        opts.max_mismatches,
        id_offset,
    );

    // Track existing 0-bulge hits by (strand, pam_start_pos) to avoid duplicate or inferior bulge calls
    let mut exact_hit_pams = std::collections::HashSet::new();
    for s in &exact_sites {
        exact_hit_pams.insert((s.strand, s.start, s.end));
    }
    sites.extend(exact_sites);

    // If no bulge search requested, return 0-bulge hits immediately
    if opts.max_dna_bulge == 0 && opts.max_rna_bulge == 0 {
        return sites;
    }

    match profile.pam_side {
        PamSide::ThreePrime => {
            // SpCas9-style: 3' PAM
            // Plus strand
            for p in 0..=seq.len().saturating_sub(pn) {
                let pam_slice = &seq[p..p + pn];
                if !iupac_match_slice(pam_bytes, pam_slice) {
                    continue;
                }
                let pam_str = String::from_utf8_lossy(pam_slice).into_owned();

                // 1-bp DNA bulge: protospacer has length sp + 1
                if opts.max_dna_bulge >= 1 && p > sp {
                    let proto_slice = &seq[p - (sp + 1)..p];
                    let mm = min_mismatches_dna_bulge_1(proto_slice, guide_bytes);
                    let start = (p - (sp + 1)) as u64 + 1;
                    let end = (p + pn) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Plus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Plus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p - (sp + 1)..p + pn])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Dna,
                            bulge_size: 1,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }

                // 2-bp DNA bulge: protospacer has length sp + 2
                if opts.max_dna_bulge >= 2 && p >= sp + 2 {
                    let proto_slice = &seq[p - (sp + 2)..p];
                    let mm = min_mismatches_dna_bulge_2(proto_slice, guide_bytes);
                    let start = (p - (sp + 2)) as u64 + 1;
                    let end = (p + pn) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Plus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Plus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p - (sp + 2)..p + pn])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Dna,
                            bulge_size: 2,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }

                // 1-bp RNA bulge: protospacer has length sp - 1
                if opts.max_rna_bulge >= 1 && p >= sp.saturating_sub(1) {
                    let proto_slice = &seq[p - (sp - 1)..p];
                    let mm = min_mismatches_rna_bulge_1(proto_slice, guide_bytes);
                    let start = (p - (sp - 1)) as u64 + 1;
                    let end = (p + pn) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Plus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Plus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p - (sp - 1)..p + pn])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Rna,
                            bulge_size: 1,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }

                // 2-bp RNA bulge: protospacer has length sp - 2
                if opts.max_rna_bulge >= 2 && p >= sp.saturating_sub(2) {
                    let proto_slice = &seq[p - (sp - 2)..p];
                    let mm = min_mismatches_rna_bulge_2(proto_slice, guide_bytes);
                    let start = (p - (sp - 2)) as u64 + 1;
                    let end = (p + pn) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Plus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Plus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p - (sp - 2)..p + pn])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Rna,
                            bulge_size: 2,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }
            }

            // Minus strand
            for p in 0..=seq.len().saturating_sub(pn) {
                let pam_slice = &seq[p..p + pn];
                if !iupac_match_slice(&pam_rc, pam_slice) {
                    continue;
                }
                let pam_rc_bytes = revcomp_dna(pam_slice);
                let pam_str = String::from_utf8_lossy(&pam_rc_bytes).into_owned();

                // 1-bp DNA bulge
                if opts.max_dna_bulge >= 1 && p + pn + sp < seq.len() {
                    let proto_rc = &seq[p + pn..p + pn + sp + 1];
                    let proto = revcomp_dna(proto_rc);
                    let mm = min_mismatches_dna_bulge_1(&proto, guide_bytes);
                    let start = p as u64 + 1;
                    let end = (p + pn + sp + 1) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Minus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Minus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p..p + pn + sp + 1])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Dna,
                            bulge_size: 1,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }

                // 1-bp RNA bulge
                if opts.max_rna_bulge >= 1 && p + pn + sp.saturating_sub(1) <= seq.len() {
                    let proto_rc = &seq[p + pn..p + pn + sp - 1];
                    let proto = revcomp_dna(proto_rc);
                    let mm = min_mismatches_rna_bulge_1(&proto, guide_bytes);
                    let start = p as u64 + 1;
                    let end = (p + pn + sp - 1) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Minus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Minus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p..p + pn + sp - 1])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Rna,
                            bulge_size: 1,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }
            }
        }
        PamSide::FivePrime => {
            // Cas12a-style: 5' PAM
            // Plus strand: [PAM][protospacer]
            for p in 0..=seq.len().saturating_sub(pn) {
                let pam_slice = &seq[p..p + pn];
                if !iupac_match_slice(pam_bytes, pam_slice) {
                    continue;
                }
                let pam_str = String::from_utf8_lossy(pam_slice).into_owned();

                // 1-bp DNA bulge
                if opts.max_dna_bulge >= 1 && p + pn + sp < seq.len() {
                    let proto_slice = &seq[p + pn..p + pn + sp + 1];
                    let mm = min_mismatches_dna_bulge_1(proto_slice, guide_bytes);
                    let start = p as u64 + 1;
                    let end = (p + pn + sp + 1) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Plus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Plus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p..p + pn + sp + 1])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Dna,
                            bulge_size: 1,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }

                // 1-bp RNA bulge
                if opts.max_rna_bulge >= 1 && p + pn + sp.saturating_sub(1) <= seq.len() {
                    let proto_slice = &seq[p + pn..p + pn + sp - 1];
                    let mm = min_mismatches_rna_bulge_1(proto_slice, guide_bytes);
                    let start = p as u64 + 1;
                    let end = (p + pn + sp - 1) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Plus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Plus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p..p + pn + sp - 1])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Rna,
                            bulge_size: 1,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }
            }

            // Minus strand: [revcomp(protospacer)][revcomp(PAM)]
            for p in 0..=seq.len().saturating_sub(pn) {
                let pam_slice = &seq[p..p + pn];
                if !iupac_match_slice(&pam_rc, pam_slice) {
                    continue;
                }
                let pam_rc_bytes = revcomp_dna(pam_slice);
                let pam_str = String::from_utf8_lossy(&pam_rc_bytes).into_owned();

                // 1-bp DNA bulge
                if opts.max_dna_bulge >= 1 && p > sp {
                    let proto_rc = &seq[p - (sp + 1)..p];
                    let proto = revcomp_dna(proto_rc);
                    let mm = min_mismatches_dna_bulge_1(&proto, guide_bytes);
                    let start = (p - (sp + 1)) as u64 + 1;
                    let end = (p + pn) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Minus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Minus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p - (sp + 1)..p + pn])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Dna,
                            bulge_size: 1,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }

                // 1-bp RNA bulge
                if opts.max_rna_bulge >= 1 && p >= sp.saturating_sub(1) {
                    let proto_rc = &seq[p - (sp - 1)..p];
                    let proto = revcomp_dna(proto_rc);
                    let mm = min_mismatches_rna_bulge_1(&proto, guide_bytes);
                    let start = (p - (sp - 1)) as u64 + 1;
                    let end = (p + pn) as u64;
                    if mm <= opts.max_mismatches
                        && !exact_hit_pams.contains(&(Strand::Minus, start, end))
                    {
                        let site_idx = id_offset + sites.len() + 1;
                        sites.push(OffTargetSite {
                            site_id: format!("site_{site_idx}"),
                            seq_id: seq_id.to_string(),
                            start,
                            end,
                            strand: Strand::Minus,
                            guide: guide_clean.clone(),
                            target_seq: String::from_utf8_lossy(&seq[p - (sp - 1)..p + pn])
                                .into_owned(),
                            pam: pam_str.clone(),
                            mismatches: mm,
                            bulge_type: BulgeType::Rna,
                            bulge_size: 1,
                            search_backend: "rust_bulge".to_string(),
                            cfd_score: None,
                            hsu_score: None,
                        });
                    }
                }
            }
        }
    }

    sites
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_min_mismatches_dna_bulge_1() {
        // Guide is 10 bp: "ACGTACGTAC"
        let guide = b"ACGTACGTAC";
        // Target has 1 inserted 'T' at index 4: "ACGT TACGTAC" (11 bp)
        let target = b"ACGTTACGTAC";
        assert_eq!(min_mismatches_dna_bulge_1(target, guide), 0);

        // Target with bulge + 1 mismatch ('G' -> 'C' at index 8)
        let target_mm = b"ACGTTACGTCC";
        assert_eq!(min_mismatches_dna_bulge_1(target_mm, guide), 1);
    }

    #[test]
    fn test_min_mismatches_rna_bulge_1() {
        // Guide is 10 bp: "ACGTACGTAC"
        let guide = b"ACGTACGTAC";
        // Target is 9 bp, missing base 4 ('A'): "ACGT CGTAC"
        let target = b"ACGTCGTAC";
        assert_eq!(min_mismatches_rna_bulge_1(target, guide), 0);

        // Target with RNA bulge + 1 mismatch
        let target_mm = b"ACGTCGTCC";
        assert_eq!(min_mismatches_rna_bulge_1(target_mm, guide), 1);
    }

    #[test]
    fn test_scan_contig_finds_dna_and_rna_bulge_sites() {
        let spacer = "GAGTCCGAGCAGAAGAAGAA"; // 20 bp
        let profile = NucleaseProfile::spcas9();

        // Target 1: 1-bp DNA bulge (21 bp protospacer + NGG)
        // Insert 'T' at position 10: "GAGTCCGAGCTAAGAAGAAGAA" (21 bp) + "CGG" (PAM)
        let mut ref_seq = vec![b'A'; 200];
        let mut proto_dna_bulge = spacer.to_string();
        proto_dna_bulge.insert(10, 'T');
        let full_dna_target = format!("{}CGG", proto_dna_bulge);
        ref_seq[20..20 + 24].copy_from_slice(full_dna_target.as_bytes());

        // Target 2: 1-bp RNA bulge (19 bp protospacer + NGG)
        // Delete base at position 5: "GAGTC GAGCAGAAGAAGAA" (19 bp) + "TGG" (PAM)
        let mut proto_rna_bulge = spacer.to_string();
        proto_rna_bulge.remove(5);
        let full_rna_target = format!("{}TGG", proto_rna_bulge);
        ref_seq[80..80 + 22].copy_from_slice(full_rna_target.as_bytes());

        let opts = BulgeSearchOptions {
            max_mismatches: 1,
            max_dna_bulge: 1,
            max_rna_bulge: 1,
        };

        let sites = scan_contig_bulge("chr1", &ref_seq, spacer, &profile, opts, 0);

        let dna_site = sites.iter().find(|s| s.bulge_type == BulgeType::Dna);
        assert!(dna_site.is_some(), "should find DNA bulge site");
        let ds = dna_site.unwrap();
        assert_eq!(ds.start, 21);
        assert_eq!(ds.end, 44);
        assert_eq!(ds.bulge_size, 1);
        assert_eq!(ds.mismatches, 0);
        assert_eq!(ds.pam, "CGG");

        let rna_site = sites.iter().find(|s| s.bulge_type == BulgeType::Rna);
        assert!(rna_site.is_some(), "should find RNA bulge site");
        let rs = rna_site.unwrap();
        assert_eq!(rs.start, 81);
        assert_eq!(rs.end, 102);
        assert_eq!(rs.bulge_size, 1);
        assert_eq!(rs.mismatches, 0);
        assert_eq!(rs.pam, "TGG");
    }
}
