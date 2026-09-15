use crate::model::{BulgeType, OffTargetSite, Strand};
use crate::normalize::clean_nucleotide_string;

#[derive(Debug, thiserror::Error)]
pub enum CasOffinderV3ParseError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error on line {line}: {msg}")]
    Format { line: usize, msg: String },
}

/// EXPERIMENTAL: Parses Cas-OFFinder 3.x extended format with bulge annotations.
pub fn parse_cas_offinder_v3(text: &str) -> Result<Vec<OffTargetSite>, CasOffinderV3ParseError> {
    let mut sites = Vec::new();
    for (line_idx, line) in text.lines().enumerate() {
        let line_num = line_idx + 1;
        let line_s = line.trim();
        if line_s.is_empty() || line_s.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line_s.split('\t').map(|s| s.trim()).collect();
        if cols.len() < 6 {
            return Err(CasOffinderV3ParseError::Format {
                line: line_num,
                msg: format!("expected at least 6 columns, got {}", cols.len()),
            });
        }

        let guide = clean_nucleotide_string(cols[0]);
        let seq_id = cols[1].to_string();
        let loc_0based: u64 = cols[2]
            .parse()
            .map_err(|_| CasOffinderV3ParseError::Format {
                line: line_num,
                msg: format!("invalid 0-based location '{}'", cols[2]),
            })?;
        let target_seq = clean_nucleotide_string(cols[3]);
        let strand = cols[4]
            .parse::<Strand>()
            .map_err(|e| CasOffinderV3ParseError::Format {
                line: line_num,
                msg: format!("invalid strand '{}': {e}", cols[4]),
            })?;
        let mismatches: u32 = cols[5]
            .parse()
            .map_err(|_| CasOffinderV3ParseError::Format {
                line: line_num,
                msg: format!("invalid mismatch count '{}'", cols[5]),
            })?;

        let mut bulge_type = BulgeType::None;
        let mut bulge_size = 0;
        if cols.len() >= 8 {
            bulge_type = cols[6].parse::<BulgeType>().unwrap_or(BulgeType::None);
            bulge_size = cols[7].parse::<u32>().unwrap_or(0);
        }

        let len = target_seq.len() as u64;
        let start = loc_0based + 1;
        let end = loc_0based + len;

        let pam = if target_seq.len() > guide.len() {
            target_seq[guide.len()..].to_string()
        } else {
            ".".to_string()
        };

        sites.push(OffTargetSite {
            site_id: format!("cas_offinder_v3_{}", sites.len() + 1),
            seq_id,
            start,
            end,
            strand,
            guide,
            target_seq,
            pam,
            mismatches,
            bulge_type,
            bulge_size,
            search_backend: "cas-offinder-v3".to_string(),
            cfd_score: None,
            hsu_score: None,
        });
    }
    Ok(sites)
}
