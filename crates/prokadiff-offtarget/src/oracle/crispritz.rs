use crate::model::{BulgeType, OffTargetSite, Strand};
use crate::normalize::clean_nucleotide_string;

#[derive(Debug, thiserror::Error)]
pub enum CrispritzParseError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error on line {line}: {msg}")]
    Format { line: usize, msg: String },
}

/// Parses CRISPRitz target output files.
pub fn parse_crispritz(text: &str) -> Result<Vec<OffTargetSite>, CrispritzParseError> {
    let mut lines = text.lines().enumerate().filter(|(_, l)| {
        let t = l.trim();
        !t.is_empty() && !t.starts_with('#')
    });

    let (hdr_idx, hdr_line) = lines.next().ok_or_else(|| CrispritzParseError::Format {
        line: 0,
        msg: "empty CRISPRitz output".into(),
    })?;

    let headers: Vec<String> = hdr_line
        .split('\t')
        .map(|s| s.trim().to_ascii_lowercase())
        .collect();

    let find_col = |names: &[&str]| -> Option<usize> {
        for name in names {
            if let Some(pos) = headers.iter().position(|h| h == name) {
                return Some(pos);
            }
        }
        None
    };

    let col_chr =
        find_col(&["chr", "chromosome", "seq_id"]).ok_or_else(|| CrispritzParseError::Format {
            line: hdr_idx + 1,
            msg: "missing chr/chromosome column".into(),
        })?;
    let col_pos =
        find_col(&["position", "start", "pos"]).ok_or_else(|| CrispritzParseError::Format {
            line: hdr_idx + 1,
            msg: "missing position/start column".into(),
        })?;
    let col_stop = find_col(&["stop", "end"]);
    let col_target = find_col(&["target", "sequence", "target_seq"]).ok_or_else(|| {
        CrispritzParseError::Format {
            line: hdr_idx + 1,
            msg: "missing target sequence column".into(),
        }
    })?;
    let col_pam = find_col(&["pam"]);
    let col_mm = find_col(&["mismatches", "mm"]);
    let col_bulge_type = find_col(&["bulge_type", "bulge"]);
    let col_bulge_size = find_col(&["bulge_size", "bulgesize"]);
    let col_strand =
        find_col(&["strand", "direction"]).ok_or_else(|| CrispritzParseError::Format {
            line: hdr_idx + 1,
            msg: "missing strand column".into(),
        })?;
    let col_cfd = find_col(&["cfd", "cfd_score"]);

    let mut sites = Vec::new();
    for (line_idx, line) in lines {
        let line_num = line_idx + 1;
        let parts: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();
        let get = |idx: usize| -> &str { parts.get(idx).copied().unwrap_or("") };

        let seq_id = get(col_chr).to_string();
        let start: u64 = get(col_pos)
            .parse()
            .map_err(|_| CrispritzParseError::Format {
                line: line_num,
                msg: format!("invalid position '{}'", get(col_pos)),
            })?;
        let target_seq = clean_nucleotide_string(get(col_target));
        let len = target_seq.len() as u64;

        let end = match col_stop {
            Some(idx) => get(idx).parse().unwrap_or(start + len.saturating_sub(1)),
            None => start + len.saturating_sub(1),
        };

        let strand =
            get(col_strand)
                .parse::<Strand>()
                .map_err(|e| CrispritzParseError::Format {
                    line: line_num,
                    msg: format!("invalid strand '{}': {e}", get(col_strand)),
                })?;

        let pam = col_pam
            .map(|i| clean_nucleotide_string(get(i)))
            .unwrap_or_else(|| ".".to_string());

        let mismatches: u32 = col_mm.and_then(|i| get(i).parse().ok()).unwrap_or(0);

        let bulge_type = col_bulge_type
            .map(|i| get(i).parse::<BulgeType>().unwrap_or(BulgeType::None))
            .unwrap_or(BulgeType::None);

        let bulge_size: u32 = col_bulge_size
            .and_then(|i| get(i).parse().ok())
            .unwrap_or(0);

        let cfd_score: Option<f64> = col_cfd.and_then(|i| {
            let s = get(i);
            if s.is_empty() || s.eq_ignore_ascii_case("na") {
                None
            } else {
                s.parse().ok()
            }
        });

        sites.push(OffTargetSite {
            site_id: format!("crispritz_{}", sites.len() + 1),
            seq_id,
            start,
            end,
            strand,
            guide: "".to_string(),
            target_seq,
            pam,
            mismatches,
            bulge_type,
            bulge_size,
            search_backend: "crispritz".to_string(),
            cfd_score,
            hsu_score: None,
        });
    }

    Ok(sites)
}
