/// Hsu et al. (2013) MIT specificity score calculation for SpCas9 (20 nt spacer).
///
/// Formula:
///   Score = ( \prod_{i \in M} (1 - W[i]) ) * ( 1 / ( ((19 - d)/19) * 4 + 1 ) ) * ( 1 / n^2 ) * 100
/// where:
///   - M is the set of mismatch indices (1-based from 5' distal to 3' PAM-proximal, 1..=20)
///   - d is the mean pairwise distance between consecutive mismatches
///   - n is the total number of mismatches
///   - if n == 0, Score = 100.0
pub fn calculate_hsu_score(guide: &str, target_protospacer: &str) -> Option<f64> {
    let g = guide.as_bytes();
    let t = target_protospacer.as_bytes();
    if g.len() != 20 || t.len() != 20 {
        return None;
    }

    // Position-dependent weights from Hsu et al. (2013), Nat Biotechnol 31(9):827-832
    // Index 0 corresponds to position 1 (5' distal), index 19 corresponds to position 20 (adjacent to PAM).
    const HSU_WEIGHTS: [f64; 20] = [
        0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851,
        0.732, 0.828, 0.615, 0.804, 0.685, 0.583,
    ];

    let mut mismatch_positions = Vec::new();
    let mut weight_product = 1.0;

    for i in 0..20 {
        if !g[i].eq_ignore_ascii_case(&t[i]) {
            mismatch_positions.push(i as f64);
            weight_product *= 1.0 - HSU_WEIGHTS[i];
        }
    }

    let n = mismatch_positions.len();
    if n == 0 {
        return Some(100.0);
    }

    // Mean distance between consecutive mismatches
    let d = if n > 1 {
        let sum_diff: f64 = mismatch_positions
            .windows(2)
            .map(|w| (w[1] - w[0]).abs())
            .sum();
        sum_diff / (n - 1) as f64
    } else {
        19.0
    };

    let dist_factor = 1.0 / (((19.0 - d) / 19.0) * 4.0 + 1.0);
    let count_factor = 1.0 / (n as f64 * n as f64);

    let score = weight_product * dist_factor * count_factor * 100.0;
    Some(score.clamp(0.0, 100.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_match_has_100_hsu_score() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        let score = calculate_hsu_score(guide, guide).unwrap();
        assert!((score - 100.0).abs() < 1e-6);
    }

    #[test]
    fn single_distal_mismatch_has_high_score() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        let mut target = guide.to_string();
        target.replace_range(0..1, "A"); // pos 1 distal mismatch (weight 0)
        let score = calculate_hsu_score(guide, &target).unwrap();
        assert!(score > 90.0);
    }

    #[test]
    fn single_seed_mismatch_has_low_score() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        let mut target = guide.to_string();
        target.replace_range(17..18, "T"); // pos 18 proximal seed mismatch (weight 0.804)
        let score = calculate_hsu_score(guide, &target).unwrap();
        assert!(score < 30.0);
    }
}
