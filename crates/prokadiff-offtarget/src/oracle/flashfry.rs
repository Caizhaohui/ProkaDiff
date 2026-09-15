use crate::model::{BulgeType, OffTargetSite, Strand};
use crate::normalize::clean_nucleotide_string;

#[derive(Debug, thiserror::Error)]
pub enum FlashFryParseError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error on line {line}: {msg}")]
    Format { line: usize, msg: String },
}

/// Parses FlashFry output TSV files with flexible header inspection.
pub fn parse_flashfry(text: &str) -> Result<Vec<OffTargetSite>, FlashFryParseError> {
    let mut lines = text.lines().enumerate().filter(|(_, l)| {
        let t = l.trim();
        !t.is_empty() && !t.starts_with('#')
    });

    let (hdr_idx, hdr_line) = lines.next().ok_or_else(|| FlashFryParseError::Format {
        line: 0,
        msg: "empty FlashFry output".into(),
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

    let col_contig = find_col(&["contig", "chrom", "chr", "seq_id"]).ok_or_else(|| {
        FlashFryParseError::Format {
            line: hdr_idx + 1,
            msg: "missing contig/chromosome column".into(),
        }
    })?;
    let col_start = find_col(&["start", "start_pos", "pos", "position"]).ok_or_else(|| {
        FlashFryParseError::Format {
            line: hdr_idx + 1,
            msg: "missing start coordinate column".into(),
        }
    })?;
    let col_stop = find_col(&["stop", "end", "end_pos"]);
    let col_strand = find_col(&["strand", "direction", "orientation"]).ok_or_else(|| {
        FlashFryParseError::Format {
            line: hdr_idx + 1,
            msg: "missing strand column".into(),
        }
    })?;
    let col_target =
        find_col(&["target", "sequence", "off_target", "target_seq"]).ok_or_else(|| {
            FlashFryParseError::Format {
                line: hdr_idx + 1,
                msg: "missing target sequence column".into(),
            }
        })?;
    let col_guide = find_col(&["guide", "crrna", "spacer"]);
    let col_mm = find_col(&["mismatches", "mismatch_count", "mm"]);
    let col_cfd = find_col(&["cfd", "cfd_score", "doench2016cfd"]);
    let col_hsu = find_col(&["hsu", "hsu_score", "hsu2013"]);

    let mut sites = Vec::new();
    for (line_idx, line) in lines {
        let line_num = line_idx + 1;
        let parts: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();
        let get = |idx: usize| -> &str { parts.get(idx).copied().unwrap_or("") };

        let seq_id = get(col_contig).to_string();
        let raw_start: u64 = get(col_start)
            .parse()
            .map_err(|_| FlashFryParseError::Format {
                line: line_num,
                msg: format!("invalid start coordinate '{}'", get(col_start)),
            })?;
        let target_seq = clean_nucleotide_string(get(col_target));
        let len = target_seq.len() as u64;

        let start = raw_start;
        let end = match col_stop {
            Some(idx) => get(idx).parse().unwrap_or(start + len.saturating_sub(1)),
            None => start + len.saturating_sub(1),
        };

        let strand = get(col_strand)
            .parse::<Strand>()
            .map_err(|e| FlashFryParseError::Format {
                line: line_num,
                msg: format!("invalid strand '{}': {e}", get(col_strand)),
            })?;

        let guide = col_guide
            .map(|i| clean_nucleotide_string(get(i)))
            .unwrap_or_default();

        let mismatches: u32 = col_mm.and_then(|i| get(i).parse().ok()).unwrap_or(0);

        let cfd_score: Option<f64> = col_cfd.and_then(|i| {
            let s = get(i);
            if s.is_empty() || s.eq_ignore_ascii_case("na") {
                None
            } else {
                s.parse().ok()
            }
        });

        let hsu_score: Option<f64> = col_hsu.and_then(|i| {
            let s = get(i);
            if s.is_empty() || s.eq_ignore_ascii_case("na") {
                None
            } else {
                s.parse().ok()
            }
        });

        let pam = if target_seq.len() > guide.len() && !guide.is_empty() {
            target_seq[guide.len()..].to_string()
        } else {
            ".".to_string()
        };

        sites.push(OffTargetSite {
            site_id: format!("flashfry_{}", sites.len() + 1),
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
            search_backend: "flashfry".to_string(),
            cfd_score,
            hsu_score,
        });
    }

    Ok(sites)
}
