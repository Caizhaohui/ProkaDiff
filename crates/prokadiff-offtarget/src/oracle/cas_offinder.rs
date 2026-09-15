pub use crate::oracle::cas_offinder_v2::{parse_cas_offinder_v2, CasOffinderV2ParseError};

/// Alias for `parse_cas_offinder_v2` for standard Cas-OFFinder mismatch-only output.
pub fn parse_cas_offinder(
    text: &str,
) -> Result<Vec<crate::model::OffTargetSite>, CasOffinderV2ParseError> {
    parse_cas_offinder_v2(text)
}
