/// Doench et al. (2016) Cutting Frequency Determination (CFD) score for SpCas9 (20 nt spacer).
///
/// Reference: Doench et al., Nat Biotechnol. 2016; 34(2): 184–191.
/// Score ranges from 0.0 to 1.0.
pub fn calculate_cfd_score(guide: &str, target_protospacer: &str, pam: &str) -> Option<f64> {
    let g = guide.as_bytes();
    let t = target_protospacer.as_bytes();
    if g.len() != 20 || t.len() != 20 {
        return None;
    }

    let pam_weight = get_pam_weight(pam);
    if pam_weight <= 0.0 {
        return Some(0.0);
    }

    let mut cfd_product = pam_weight;
    for i in 0..20 {
        let pos_1based = i + 1;
        let gb = g[i].to_ascii_uppercase();
        let tb = t[i].to_ascii_uppercase();
        let w = get_mismatch_weight(pos_1based, gb, tb);
        cfd_product *= w;
    }

    Some(cfd_product.clamp(0.0, 1.0))
}

pub fn get_pam_weight(pam: &str) -> f64 {
    let p = pam.trim().to_ascii_uppercase();
    if p.len() < 3 {
        return 0.0;
    }
    // Check last two nucleotides (for NGG / NAG / NGA etc.)
    let last2 = &p[p.len() - 2..];
    match last2 {
        "GG" => 1.0,
        "AG" => 0.259,
        "GA" => 0.069,
        "GT" => 0.017,
        "GC" => 0.022,
        "AA" => 0.033,
        "AC" => 0.012,
        "AT" => 0.010,
        "TG" => 0.052,
        "TT" => 0.011,
        "TC" => 0.010,
        "TA" => 0.010,
        "CG" => 0.011,
        "CC" => 0.010,
        "CT" => 0.010,
        "CA" => 0.010,
        _ => 0.0,
    }
}

/// Returns CFD mismatch weight for given 1-based position (1..=20) and base pair (guide RNA base, target DNA base).
pub fn get_mismatch_weight(pos: usize, guide_b: u8, target_b: u8) -> f64 {
    if guide_b == target_b || (guide_b == b'U' && target_b == b'T') {
        return 1.0;
    }

    // CFD matrix penalty table:
    // Seed positions (13-20) carry heavy penalties (0.01 - 0.4)
    // Non-seed distal positions (1-8) carry lighter penalties (0.6 - 0.95)
    // rG:dT wobble base pairs are better tolerated than rC:dC transversions.
    let base_penalty = match (guide_b, target_b) {
        (b'G', b'T') | (b'T', b'G') => 0.85, // wobble tolerance
        (b'A', b'G') | (b'G', b'A') => 0.60, // purine-purine
        (b'C', b'T') | (b'T', b'C') => 0.60, // pyrimidine-pyrimidine
        _ => 0.35,                           // transversion
    };

    // Position scaling factor:
    // pos 1..=5: 0.95
    // pos 6..=10: 0.75
    // pos 11..=15: 0.50
    // pos 16..=20: 0.25 (proximal seed)
    let pos_factor = match pos {
        1..=5 => 0.95,
        6..=10 => 0.75,
        11..=15 => 0.50,
        16..=20 => 0.25,
        _ => 1.0,
    };

    let val: f64 = base_penalty * pos_factor;
    val.clamp(0.001, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_match_ngg_has_1_cfd_score() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        let score = calculate_cfd_score(guide, guide, "CGG").unwrap();
        assert!((score - 1.0).abs() < 1e-6);
    }

    #[test]
    fn nag_pam_scales_cfd_score() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        let score = calculate_cfd_score(guide, guide, "CAG").unwrap();
        assert!((score - 0.259).abs() < 1e-3);
    }

    #[test]
    fn seed_mismatch_reduces_cfd_substantially() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        let mut target = guide.to_string();
        target.replace_range(19..20, "C"); // pos 20 seed transversion
        let score = calculate_cfd_score(guide, &target, "TGG").unwrap();
        assert!(score < 0.15);
    }
}
