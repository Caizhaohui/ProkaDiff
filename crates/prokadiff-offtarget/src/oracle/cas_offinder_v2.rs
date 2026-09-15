use crate::model::{BulgeType, OffTargetSite, Strand};
use crate::normalize::clean_nucleotide_string;

#[derive(Debug, thiserror::Error)]
pub enum CasOffinderV2ParseError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error on line {line}: {msg}")]
    Format { line: usize, msg: String },
}

/// Parses Cas-OFFinder 2.x standard 6-column tab-delimited mismatch-only output.
///
/// Schema:
///   `CrRNA \t Chromosome \t Location(0-based) \t Target_Sequence \t Direction(+/-) \t Mismatches`
pub fn parse_cas_offinder_v2(text: &str) -> Result<Vec<OffTargetSite>, CasOffinderV2ParseError> {
    let mut sites = Vec::new();
    for (line_idx, line) in text.lines().enumerate() {
        let line_num = line_idx + 1;
        let line_s = line.trim();
        if line_s.is_empty() || line_s.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line_s.split('\t').map(|s| s.trim()).collect();
        if cols.len() < 6 {
            return Err(CasOffinderV2ParseError::Format {
                line: line_num,
                msg: format!(
                    "expected exactly 6 columns in Cas-OFFinder v2 output, got {}",
                    cols.len()
                ),
            });
        }

        let guide = clean_nucleotide_string(cols[0]);
        let seq_id = cols[1].to_string();
        let loc_0based: u64 = cols[2]
            .parse()
            .map_err(|_| CasOffinderV2ParseError::Format {
                line: line_num,
                msg: format!("invalid 0-based location '{}'", cols[2]),
            })?;
        let target_seq = clean_nucleotide_string(cols[3]);
        let strand = cols[4]
            .parse::<Strand>()
            .map_err(|e| CasOffinderV2ParseError::Format {
                line: line_num,
                msg: format!("invalid strand '{}': {e}", cols[4]),
            })?;
        let mismatches: u32 = cols[5]
            .parse()
            .map_err(|_| CasOffinderV2ParseError::Format {
                line: line_num,
                msg: format!("invalid mismatch count '{}'", cols[5]),
            })?;

        let len = target_seq.len() as u64;
        let start = loc_0based + 1; // 1-based start
        let end = loc_0based + len; // 1-based inclusive end

        let pam = if target_seq.len() > guide.len() {
            target_seq[guide.len()..].to_string()
        } else {
            ".".to_string()
        };

        sites.push(OffTargetSite {
            site_id: format!("cas_offinder_{}", sites.len() + 1),
            seq_id,
            start,
            end,
            strand,
            guide,
            target_seq,
            pam,
            mismatches,
            bulge_type: BulgeType::None,
            bulge_size: 0,
            search_backend: "cas-offinder-v2".to_string(),
            cfd_score: None,
            hsu_score: None,
        });
    }
    Ok(sites)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parser_accepts_known_v2_fixture() {
        let sample = "\
GAGTCCGAGCAGAAGAAGAA\tchr1\t100\tGAGTCCGAGCAGAAGAAGAATGG\t+\t0
GAGTCCGAGCAGAAGAAGAA\tchr1\t250\tGAGTCCGAGCAGAAGAACAAGGG\t-\t1
";
        let sites = parse_cas_offinder_v2(sample).expect("should parse Cas-OFFinder v2 output");
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].start, 101);
        assert_eq!(sites[0].end, 123);
        assert_eq!(sites[0].strand, Strand::Plus);
        assert_eq!(sites[0].mismatches, 0);
        assert_eq!(sites[0].bulge_type, BulgeType::None);
        assert_eq!(sites[0].bulge_size, 0);
        assert_eq!(sites[1].strand, Strand::Minus);
        assert_eq!(sites[1].mismatches, 1);
    }
}
