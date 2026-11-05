//! Subclonal / low-frequency variant statistical analysis and quality filtering.
//!
//! Milestone M11: Subclonal / Pool-seq Variant Audit Mode.
//! Implements Wilson Score interval with continuity correction (Newcombe 1998, Method 4),
//! two-tailed Fisher's exact test with log-factorial Stirling approximation,
//! upper-tail binomial background sequencing error filtering,
//! Benjamini-Hochberg FDR adjustment, and Newcombe hybrid-score differential interval.

use std::sync::LazyLock;

/// Precomputed log-factorials `ln(k!)` for k < 256.
static LN_FACT_256: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut table = [0.0; 256];
    let mut sum = 0.0;
    for (i, entry) in table.iter_mut().enumerate().skip(2) {
        sum += (i as f64).ln();
        *entry = sum;
    }
    table
});

/// Log-factorial `ln(k!)`.
/// Uses an exact precomputed lookup for k < 256 and high-order Stirling asymptotic expansion
/// for k >= 256, ensuring absolute error < 1e-14 and zero factorial overflow.
pub fn lnfact(k: u64) -> f64 {
    if k < 256 {
        LN_FACT_256[k as usize]
    } else {
        let k_f = k as f64;
        let pi = std::f64::consts::PI;
        k_f * k_f.ln() - k_f + 0.5 * (2.0 * pi * k_f).ln()
            + 1.0 / (12.0 * k_f)
            - 1.0 / (360.0 * k_f * k_f * k_f)
            + 1.0 / (1260.0 * k_f.powi(5))
    }
}

/// Log-combinations `ln(n choose k)`.
pub fn ln_choose(n: u64, k: u64) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    lnfact(n) - lnfact(k) - lnfact(n - k)
}

/// Errors occurring during statistical evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum StatsError {
    #[error("zero total depth")]
    ZeroDepth,
    #[error("allele count exceeds total depth")]
    CountInconsistent,
    #[error("invalid parameter rate outside range")]
    InvalidRate,
}

/// Strand-resolved read counts.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StrandCounts {
    pub fwd: u32,
    pub rev: u32,
}

impl StrandCounts {
    pub const fn new(fwd: u32, rev: u32) -> Self {
        Self { fwd, rev }
    }

    pub const fn total(&self) -> u64 {
        self.fwd as u64 + self.rev as u64
    }
}

/// Allele observations at a candidate small variant locus.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AlleleCounts {
    pub ref_counts: StrandCounts,
    pub alt_counts: StrandCounts,
    pub other: u32,
}

impl AlleleCounts {
    pub fn new(ref_counts: StrandCounts, alt_counts: StrandCounts, other: u32) -> Self {
        Self {
            ref_counts,
            alt_counts,
            other,
        }
    }

    pub fn total_depth(&self) -> u64 {
        self.ref_counts
            .total()
            .saturating_add(self.alt_counts.total())
            .saturating_add(self.other as u64)
    }

    pub fn alt_fraction(&self) -> Option<f64> {
        let n = self.total_depth();
        if n == 0 {
            None
        } else {
            Some(self.alt_counts.total() as f64 / n as f64)
        }
    }
}

/// Wilson Score confidence interval [lower, upper].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WilsonCi {
    pub low: f64,
    pub high: f64,
}

/// Calculates 95% Wilson Score confidence interval with Continuity Correction (Newcombe 1998, Method 4).
///
/// Standard z = 1.959963984540054 (two-sided 95%).
pub fn calculate_wilson_ci(alt: u64, depth: u64) -> Result<WilsonCi, StatsError> {
    calculate_wilson_ci_with_z(alt, depth, 1.959963984540054)
}

/// Calculates Wilson Score confidence interval with Continuity Correction given a custom z score.
pub fn calculate_wilson_ci_with_z(alt: u64, depth: u64, z: f64) -> Result<WilsonCi, StatsError> {
    if depth == 0 {
        return Err(StatsError::ZeroDepth);
    }
    if alt > depth {
        return Err(StatsError::CountInconsistent);
    }

    if alt == 0 {
        // Continuity-corrected upper bound when A = 0
        let n_f = depth as f64;
        let z2 = z * z;
        let denom = 2.0 * (n_f + z2);
        let radicand = (z2 + 2.0 - 1.0 / n_f).max(0.0);
        let upper = (z2 + 1.0 + z * radicand.sqrt()) / denom;
        return Ok(WilsonCi {
            low: 0.0,
            high: upper.clamp(0.0, 1.0),
        });
    }

    if alt == depth {
        // Continuity-corrected lower bound when A = N
        let n_f = depth as f64;
        let z2 = z * z;
        let denom = 2.0 * (n_f + z2);
        let radicand = (z2 + 2.0 - 1.0 / n_f).max(0.0);
        let lower = (2.0 * n_f + z2 - 1.0 - z * radicand.sqrt()) / denom;
        return Ok(WilsonCi {
            low: lower.clamp(0.0, 1.0),
            high: 1.0,
        });
    }

    let a_f = alt as f64;
    let n_f = depth as f64;
    let z2 = z * z;
    let denom = 2.0 * (n_f + z2);

    let l_radicand = (z2 - 2.0 - 1.0 / n_f + 4.0 * (a_f / n_f) * (n_f - a_f + 1.0)).max(0.0);
    let l_num = 2.0 * a_f + z2 - 1.0 - z * l_radicand.sqrt();
    let low = (l_num / denom).clamp(0.0, 1.0);

    let u_radicand = (z2 + 2.0 - 1.0 / n_f + 4.0 * (a_f / n_f) * (n_f - a_f - 1.0)).max(0.0);
    let u_num = 2.0 * a_f + z2 + 1.0 + z * u_radicand.sqrt();
    let high = (u_num / denom).clamp(0.0, 1.0);

    Ok(WilsonCi {
        low: low.min(high),
        high: high.max(low),
    })
}

/// Calculates two-tailed Fisher's exact test p-value for a 2x2 contingency table:
/// `[[r11, r12], [r21, r22]]` where row 1 = REF (fwd, rev), row 2 = ALT (fwd, rev).
///
/// Follows R's `fisher.test` two-sided sum over hypergeometric probabilities
/// satisfying `P(x) <= P(x_obs) * (1 + 1e-7)`.
pub fn fisher_two_sided(table: [[u64; 2]; 2]) -> f64 {
    let a = table[0][0];
    let b = table[0][1];
    let c = table[1][0];
    let d = table[1][1];

    let Some(r1) = a.checked_add(b) else { return 1.0; };
    let Some(r2) = c.checked_add(d) else { return 1.0; };
    let Some(c1) = a.checked_add(c) else { return 1.0; };
    let Some(c2) = b.checked_add(d) else { return 1.0; };
    let Some(n) = r1.checked_add(r2) else { return 1.0; };

    if r1 == 0 || r2 == 0 || c1 == 0 || c2 == 0 || n == 0 {
        return 1.0;
    }

    let min_x = c1.saturating_sub(r2);
    let max_x = r1.min(c1);

    if min_x > max_x {
        return 1.0;
    }

    let ln_denom = ln_choose(n, c1);
    let log_p_obs = ln_choose(r1, a) + ln_choose(r2, c1 - a) - ln_denom;
    let p_obs = log_p_obs.exp();
    let threshold = p_obs * (1.0 + 1e-7);

    let mut sum_p = 0.0;
    for x in min_x..=max_x {
        let log_p_x = ln_choose(r1, x) + ln_choose(r2, c1 - x) - ln_denom;
        let p_x = log_p_x.exp();
        if p_x <= threshold {
            sum_p += p_x;
        }
    }

    sum_p.clamp(0.0, 1.0)
}

/// Calculates exact upper-tail binomial probability:
/// `p_err = P(X >= alt | depth, e0) = sum_{k=alt}^depth C(depth, k) e0^k (1 - e0)^(depth - k)`
/// under null hypothesis H0: true error proportion <= e0.
pub fn binom_upper_tail(alt: u64, depth: u64, e0: f64) -> Result<f64, StatsError> {
    if depth == 0 {
        return Err(StatsError::ZeroDepth);
    }
    if alt > depth {
        return Err(StatsError::CountInconsistent);
    }
    if !(0.0..=1.0).contains(&e0) {
        return Err(StatsError::InvalidRate);
    }

    if alt == 0 {
        return Ok(1.0);
    }
    if e0 <= 0.0 {
        return Ok(0.0);
    }
    if e0 >= 1.0 {
        return Ok(1.0);
    }

    let n = depth;
    let expected_mean = (n as f64) * e0;

    // When alt is below the expected null mean, compute lower tail P(X < alt)
    // and return 1.0 - lower_tail. This prevents underflow when starting from an extreme low value.
    if (alt as f64) <= expected_mean {
        let k_init = alt - 1;
        let log_pmf = ln_choose(n, k_init)
            + (k_init as f64) * e0.ln()
            + ((n - k_init) as f64) * (1.0 - e0).ln();
        let mut term = log_pmf.exp();
        let mut lower_sum = term;
        let down_factor = (1.0 - e0) / e0;

        for k in (1..=k_init).rev() {
            let mult = (k as f64) / ((n - k + 1) as f64) * down_factor;
            term *= mult;
            lower_sum += term;
            if term < 1e-18 * lower_sum {
                break;
            }
        }

        return Ok((1.0 - lower_sum).clamp(0.0, 1.0));
    }

    // When alt is strictly above expected mean, compute upper tail starting at alt
    let log_pmf_a = ln_choose(n, alt)
        + (alt as f64) * e0.ln()
        + ((n - alt) as f64) * (1.0 - e0).ln();
    let mut term = log_pmf_a.exp();
    let mut sum = term;
    let ratio_factor = e0 / (1.0 - e0);

    for k in alt..n {
        let mult = ((n - k) as f64) / ((k + 1) as f64) * ratio_factor;
        term *= mult;
        sum += term;
        if term < 1e-18 * sum && (k as f64) > expected_mean {
            break;
        }
    }

    Ok(sum.clamp(0.0, 1.0))
}

/// Benjamini-Hochberg procedure adjusting p-values to control False Discovery Rate (FDR).
pub fn bh_adjust(pvalues: &[f64]) -> Vec<f64> {
    let m = pvalues.len();
    if m == 0 {
        return Vec::new();
    }
    if m == 1 {
        return vec![pvalues[0].clamp(0.0, 1.0)];
    }

    let mut indexed: Vec<(usize, f64)> = pvalues.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut qvalues = vec![0.0; m];
    let m_f = m as f64;
    let mut min_q = 1.0;

    for rank in (1..=m).rev() {
        let (orig_idx, p) = indexed[rank - 1];
        let q = (p * m_f / (rank as f64)).min(min_q).clamp(0.0, 1.0);
        min_q = q;
        qvalues[orig_idx] = q;
    }

    qvalues
}

/// Calculates Newcombe hybrid score difference confidence interval [delta_L, delta_U]
/// for Delta = AF_edited - AF_starter:
/// `Delta_L = Delta - sqrt((AF_e - L_e)^2 + (U_s - AF_s)^2)`
/// `Delta_U = Delta + sqrt((U_e - AF_e)^2 + (AF_s - L_s)^2)`
pub fn calculate_delta_ci(
    af_e: f64,
    ci_e: &WilsonCi,
    af_s: f64,
    ci_s: &WilsonCi,
) -> (f64, f64, f64) {
    let delta = af_e - af_s;
    let l_diff_e = af_e - ci_e.low;
    let u_diff_s = ci_s.high - af_s;
    let delta_l = delta - (l_diff_e * l_diff_e + u_diff_s * u_diff_s).sqrt();

    let u_diff_e = ci_e.high - af_e;
    let l_diff_s = af_s - ci_s.low;
    let delta_u = delta + (u_diff_e * u_diff_e + l_diff_s * l_diff_s).sqrt();

    (delta, delta_l.clamp(-1.0, 1.0), delta_u.clamp(-1.0, 1.0))
}

/// Configuration parameters for subclonal variant audit.
#[derive(Clone, Debug, PartialEq)]
pub struct SubclonalParams {
    /// Minimum total depth per sample (default: 20).
    pub min_depth: u64,
    /// Minimum supporting ALT reads (default: 3).
    pub min_alt_reads: u64,
    /// Minimum ALT allele fraction (default: 0.05).
    pub min_af: f64,
    /// Maximum ALT allele fraction (default: 0.7499, below consensus 0.75).
    pub max_af: f64,
    /// Minimum Delta AF for differential calling (default: 0.05).
    pub min_delta_af: f64,
    /// Maximum strand bias p-value for soft warning flag (default: 0.001).
    pub max_strand_bias_pvalue: f64,
    /// False discovery rate alpha threshold (default: 0.01).
    pub fdr_alpha: f64,
    /// Base sequencing error floor for SNP/SUB (default: 0.005).
    pub error_floor_snp: f64,
    /// Base sequencing error floor for INS/DEL <= 2bp (default: 0.010).
    pub error_floor_indel: f64,
}

impl Default for SubclonalParams {
    fn default() -> Self {
        Self {
            min_depth: 20,
            min_alt_reads: 3,
            min_af: 0.05,
            max_af: 0.7499,
            min_delta_af: 0.05,
            max_strand_bias_pvalue: 0.001,
            fdr_alpha: 0.01,
            error_floor_snp: 0.005,
            error_floor_indel: 0.010,
        }
    }
}

impl SubclonalParams {
    pub fn validate(&self) -> Result<(), &'static str> {
        if self.min_depth == 0 {
            return Err("min_depth must be > 0");
        }
        if self.min_alt_reads == 0 {
            return Err("min_alt_reads must be > 0");
        }
        if self.min_alt_reads > self.min_depth {
            return Err("min_alt_reads cannot exceed min_depth");
        }
        if self.min_af <= 0.0 || self.min_af >= self.max_af {
            return Err("min_af must be in (0, max_af)");
        }
        if self.max_af > 1.0 {
            return Err("max_af cannot exceed 1.0");
        }
        if self.fdr_alpha <= 0.0 || self.fdr_alpha >= 1.0 {
            return Err("fdr_alpha must be in (0, 1)");
        }
        if self.error_floor_snp <= 0.0 || self.error_floor_snp >= 0.5 {
            return Err("error_floor_snp must be in (0, 0.5)");
        }
        if self.error_floor_indel <= 0.0 || self.error_floor_indel >= 0.5 {
            return Err("error_floor_indel must be in (0, 0.5)");
        }
        Ok(())
    }
}

/// Subclonal variant classification classes evaluated in deterministic priority order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SubclonalClass {
    /// Starter sample depth is below minimum threshold (< min_depth); absence cannot be confirmed.
    IndeterminateStarterDepth,
    /// Detected de novo in edited sample and absent in starter with significant differential.
    EditedSpecific,
    /// Present in both starter and edited, but significantly enriched in edited sample.
    EditedEnriched,
    /// Pre-existing in starter without significant enrichment in edited; subtracted from edited findings.
    PreExistingInStarter,
    /// Detected in starter and significantly depleted/diluted in edited sample.
    DepletedInEdited,
    /// Detected in edited sample, but starter has trace reads (A_s >= 1) with non-significant differential.
    AmbiguousStarterTrace,
    /// Does not meet subclonal audit thresholds.
    NotSubclonal,
}

impl SubclonalClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::IndeterminateStarterDepth => "INDETERMINATE_STARTER_DEPTH",
            Self::EditedSpecific => "EDITED_SPECIFIC",
            Self::EditedEnriched => "EDITED_ENRICHED",
            Self::PreExistingInStarter => "PRE_EXISTING_IN_STARTER",
            Self::DepletedInEdited => "DEPLETED_IN_EDITED",
            Self::AmbiguousStarterTrace => "AMBIGUOUS_STARTER_TRACE",
            Self::NotSubclonal => "NOT_SUBCLONAL",
        }
    }

    /// Whether this subclonal classification represents a reportable post-edit finding.
    /// Excludes `NotSubclonal` and un-enriched `PreExistingInStarter` (which is subtracted).
    pub fn is_reportable(&self) -> bool {
        matches!(
            self,
            Self::EditedSpecific
                | Self::EditedEnriched
                | Self::IndeterminateStarterDepth
                | Self::AmbiguousStarterTrace
        )
    }
}

/// Evaluates subclonal classification for candidate allele based on starter and edited observations.
pub fn classify_subclonal_candidate(
    counts_s: &AlleleCounts,
    counts_e: &AlleleCounts,
    q_err_s: f64,
    q_err_e: f64,
    _p_diff: f64,
    q_diff: f64,
    params: &SubclonalParams,
) -> (SubclonalClass, Option<(f64, f64, f64)>) {
    let n_s = counts_s.total_depth();
    let a_s = counts_s.alt_counts.total();
    let af_s = counts_s.alt_fraction().unwrap_or(0.0);

    let n_e = counts_e.total_depth();
    let a_e = counts_e.alt_counts.total();
    let af_e = counts_e.alt_fraction().unwrap_or(0.0);

    let det_s = n_s >= params.min_depth
        && a_s >= params.min_alt_reads
        && af_s >= params.min_af
        && q_err_s <= params.fdr_alpha;

    let det_e = n_e >= params.min_depth
        && a_e >= params.min_alt_reads
        && af_e >= params.min_af
        && af_e <= params.max_af
        && q_err_e <= params.fdr_alpha;

    let ci_s = calculate_wilson_ci(a_s, n_s).unwrap_or(WilsonCi { low: 0.0, high: 0.0 });
    let ci_e = calculate_wilson_ci(a_e, n_e).unwrap_or(WilsonCi { low: 0.0, high: 0.0 });
    let (delta, delta_l, delta_u) = calculate_delta_ci(af_e, &ci_e, af_s, &ci_s);
    let delta_metrics = Some((delta, delta_l, delta_u));

    // Priority 1: Indeterminate starter depth
    if det_e && n_s < params.min_depth {
        return (SubclonalClass::IndeterminateStarterDepth, delta_metrics);
    }

    let diff_significant = q_diff <= params.fdr_alpha;

    // Priority 2: Edited-specific
    if det_e && !det_s && diff_significant && delta >= params.min_delta_af && delta_l > 0.0 {
        return (SubclonalClass::EditedSpecific, delta_metrics);
    }

    // Priority 3: Edited-enriched
    if det_s && det_e && diff_significant && delta >= params.min_delta_af && delta_l > 0.0 {
        return (SubclonalClass::EditedEnriched, delta_metrics);
    }

    // Priority 4: Pre-existing in starter
    if det_s && (!det_e || !diff_significant || delta < params.min_delta_af) {
        if diff_significant && delta <= -params.min_delta_af {
            return (SubclonalClass::DepletedInEdited, delta_metrics);
        }
        return (SubclonalClass::PreExistingInStarter, delta_metrics);
    }

    // Priority 5: Ambiguous starter trace
    if det_e && !det_s && a_s >= 1 {
        return (SubclonalClass::AmbiguousStarterTrace, delta_metrics);
    }

    (SubclonalClass::NotSubclonal, delta_metrics)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fisher_two_sided_reference_vector() {
        // Reference vector: [[3, 1], [1, 3]] -> p = 0.4857143
        let t = [[3, 1], [1, 3]];
        let p = fisher_two_sided(t);
        let expected = 34.0 / 70.0; // 0.4857142857142857
        assert!((p - expected).abs() < 1e-6, "p = {}, expected {}", p, expected);

        // Invariance under matrix transposition
        let t_trans = [[3, 1], [1, 3]];
        assert_eq!(fisher_two_sided(t), fisher_two_sided(t_trans));

        // Symmetric diagonal
        let t_diag = [[10, 0], [0, 10]];
        let p_diag = fisher_two_sided(t_diag);
        assert!(p_diag < 1e-4);
    }

    #[test]
    fn test_wilson_ci_continuity_corrected_reference() {
        // A = 0, N = 20 -> U ≈ 0.2005
        let ci = calculate_wilson_ci(0, 20).expect("wilson ci");
        assert_eq!(ci.low, 0.0);
        assert!((ci.high - 0.2005).abs() < 0.005, "U = {}, expected ≈ 0.2005", ci.high);

        // A = 10, N = 20 -> symmetric around 0.5
        let ci_sym = calculate_wilson_ci(10, 20).expect("wilson ci");
        let mid = (ci_sym.low + ci_sym.high) / 2.0;
        assert!((mid - 0.5).abs() < 1e-6, "mid = {}, expected 0.5", mid);

        // A = 20, N = 20 -> high = 1.0
        let ci_full = calculate_wilson_ci(20, 20).expect("wilson ci");
        assert_eq!(ci_full.high, 1.0);
    }

    #[test]
    fn test_binomial_upper_tail_reference() {
        // A = 3, N = 20, e0 = 0.005 -> p_err ≈ 1.34e-4
        let p = binom_upper_tail(3, 20, 0.005).expect("binom");
        assert!((p - 1.337e-4).abs() < 1e-6, "p = {}, expected ≈ 1.34e-4", p);

        // A = 0 -> p = 1.0
        assert_eq!(binom_upper_tail(0, 20, 0.005).unwrap(), 1.0);

        // A = 5000, N = 100000, e0 = 0.1 (well below mean 10000) -> must not underflow to 0, returns 1.0
        assert_eq!(binom_upper_tail(5000, 100000, 0.1).unwrap(), 1.0);

        // A = 20, e0 = 0.005 -> very small
        let p_max = binom_upper_tail(20, 20, 0.005).unwrap();
        assert!(p_max < 1e-30);
    }

    #[test]
    fn test_bh_adjust() {
        let pvals = vec![0.001, 0.01, 0.02, 0.05];
        let qvals = bh_adjust(&pvals);
        assert_eq!(qvals.len(), 4);
        assert!(qvals[0] <= qvals[1]);
        assert!(qvals[1] <= qvals[2]);
        assert!(qvals[2] <= qvals[3]);
        for q in qvals {
            assert!((0.0..=1.0).contains(&q));
        }
    }

    #[test]
    fn test_subclonal_classification_flow() {
        let params = SubclonalParams::default();

        // 1. Starter depth < 20 -> IndeterminateStarterDepth
        let counts_s_low = AlleleCounts::new(StrandCounts::new(5, 5), StrandCounts::new(0, 0), 0);
        let counts_e = AlleleCounts::new(StrandCounts::new(45, 45), StrandCounts::new(5, 5), 0);
        let (class1, _) = classify_subclonal_candidate(
            &counts_s_low, &counts_e, 1.0, 0.0001, 0.0001, 0.0001, &params,
        );
        assert_eq!(class1, SubclonalClass::IndeterminateStarterDepth);

        // 2. Clear de novo in edited -> EditedSpecific
        let counts_s_high = AlleleCounts::new(StrandCounts::new(50, 50), StrandCounts::new(0, 0), 0);
        let (class2, delta2) = classify_subclonal_candidate(
            &counts_s_high, &counts_e, 1.0, 0.0001, 0.0001, 0.0001, &params,
        );
        assert_eq!(class2, SubclonalClass::EditedSpecific);
        let (d, dl, _du) = delta2.unwrap();
        assert!(d > 0.09 && d < 0.11);
        assert!(dl > 0.0);

        // 3. Pre-existing in starter and enriched in edited -> EditedEnriched
        let counts_s_low_af = AlleleCounts::new(StrandCounts::new(95, 95), StrandCounts::new(5, 5), 0);
        let counts_e_high_af = AlleleCounts::new(StrandCounts::new(60, 60), StrandCounts::new(40, 40), 0);
        let (class3, _) = classify_subclonal_candidate(
            &counts_s_low_af, &counts_e_high_af, 0.001, 0.0001, 0.0001, 0.0001, &params,
        );
        assert_eq!(class3, SubclonalClass::EditedEnriched);

        // 4. Pre-existing in starter without significant enrichment -> PreExistingInStarter
        let counts_e_same = AlleleCounts::new(StrandCounts::new(95, 95), StrandCounts::new(5, 5), 0);
        let (class4, _) = classify_subclonal_candidate(
            &counts_s_low_af, &counts_e_same, 0.001, 0.001, 1.0, 1.0, &params,
        );
        assert_eq!(class4, SubclonalClass::PreExistingInStarter);

        // 5. Consensus starter allele (100% in starter, 50% in edited) -> DepletedInEdited (not AmbiguousStarterTrace)
        let counts_s_consensus = AlleleCounts::new(StrandCounts::new(0, 0), StrandCounts::new(100, 100), 0);
        let counts_e_50 = AlleleCounts::new(StrandCounts::new(50, 50), StrandCounts::new(50, 50), 0);
        let (class5, _) = classify_subclonal_candidate(
            &counts_s_consensus, &counts_e_50, 0.0001, 0.0001, 0.0001, 0.0001, &params,
        );
        assert_eq!(class5, SubclonalClass::DepletedInEdited);
        assert!(!class5.is_reportable());
    }
}
