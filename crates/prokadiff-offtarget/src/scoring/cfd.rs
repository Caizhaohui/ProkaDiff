/// Cutting Frequency Determination (CFD) score for SpCas9 (Doench et al., 2016).
///
/// Status: Disabled pending validated FlashFry parity against official Doench 2016 matrix.
pub fn calculate_cfd_score(_guide: &str, _target_protospacer: &str, _pam: &str) -> Option<f64> {
    // Disabled pending validated FlashFry parity.
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cfd_disabled_until_validated() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        assert_eq!(calculate_cfd_score(guide, guide, "CGG"), None);
    }
}
