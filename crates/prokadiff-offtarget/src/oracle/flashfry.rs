use crate::model::{BulgeType, OffTargetSite, Strand};
use crate::normalize::clean_nucleotide_string;

/// FlashFry 1.15 native output columns (lowercase for matching):
///   contig  start  stop  target  orientation  numberOfMismatches
///   Doench2016CFDScore  Hsu2013  [other score columns...]
///
/// Coordinate convention (to be confirmed against real FlashFry 1.15 output):
///   FlashFry start: 0-based inclusive
///   FlashFry stop:  0-based exclusive (half-open)
///   → ProkaDiff: start = ff_start + 1, end = ff_stop (1-based inclusive)
///
/// Orientation:
///   FlashFry "FWD" → Strand::Plus
///   FlashFry "RVS" → Strand::Minus
///
/// Status: BLOCKED_EXTERNAL_DEPENDENCY — coordinate convention confirmed only
/// when real FlashFry 1.15 output is available.  Parser infrastructure is
/// ready; coordinate conversion follows the documented FlashFry convention.
#[derive(Debug, thiserror::Error)]
pub enum FlashFryParseError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error on line {line}: {msg}")]
    Format { line: usize, msg: String },
}

/// Parses FlashFry 1.15 output TSV.
///
/// Prioritizes native 1.15 column names; falls back to common aliases so the
/// parser tolerates minor header variations without breaking.
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

    // Resolve a column by trying a list of names in priority order.
    let find_col = |names: &[&str]| -> Option<usize> {
        for name in names {
            if let Some(pos) = headers.iter().position(|h| h == name) {
                return Some(pos);
            }
        }
        None
    };

    // ── required columns ────────────────────────────────────────────────────
    // FlashFry 1.15 native name first, then common aliases.
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
    // FlashFry 1.15 uses "orientation"; also accept "strand", "direction".
    let col_strand = find_col(&["orientation", "strand", "direction"]).ok_or_else(|| {
        FlashFryParseError::Format {
            line: hdr_idx + 1,
            msg: "missing orientation/strand column".into(),
        }
    })?;
    let col_target =
        find_col(&["target", "sequence", "off_target", "target_seq"]).ok_or_else(|| {
            FlashFryParseError::Format {
                line: hdr_idx + 1,
                msg: "missing target sequence column".into(),
            }
        })?;

    // ── optional columns ────────────────────────────────────────────────────
    let col_guide = find_col(&["guide", "crrna", "spacer"]);
    // FlashFry 1.15 mismatch column: "numberOfMismatches"
    let col_mm = find_col(&["numberofmismatches", "mismatches", "mismatch_count", "mm"]);
    // FlashFry 1.15 score columns: "Doench2016CFDScore", "Hsu2013"
    let col_cfd = find_col(&["doench2016cfdscore", "cfd", "cfd_score", "doench2016cfd"]);
    let col_hsu = find_col(&["hsu2013", "hsu", "hsu_score"]);

    let mut sites = Vec::new();
    for (line_idx, line) in lines {
        let line_num = line_idx + 1;
        let parts: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();
        let get = |idx: usize| -> &str { parts.get(idx).copied().unwrap_or("") };

        let seq_id = get(col_contig).to_string();

        // Coordinate conversion: FlashFry 0-based half-open → 1-based inclusive.
        // prokadiff_start = ff_start + 1
        // prokadiff_end   = ff_stop  (half-open exclusive → last nt at ff_stop - 1; +1 for 1-based → ff_stop)
        let ff_start: u64 = get(col_start)
            .parse()
            .map_err(|_| FlashFryParseError::Format {
                line: line_num,
                msg: format!("invalid start coordinate '{}'", get(col_start)),
            })?;

        let target_seq = clean_nucleotide_string(get(col_target));
        let len = target_seq.len() as u64;

        let start = ff_start + 1; // 0-based → 1-based
        let end = match col_stop {
            Some(idx) => {
                let ff_stop: u64 = get(idx).parse().unwrap_or(ff_start + len);
                ff_stop // half-open exclusive → 1-based inclusive end = ff_stop
            }
            None => start + len.saturating_sub(1),
        };

        let strand = get(col_strand)
            .parse::<Strand>()
            .map_err(|e| FlashFryParseError::Format {
                line: line_num,
                msg: format!("invalid orientation '{}': {e}", get(col_strand)),
            })?;

        let guide = col_guide
            .map(|i| clean_nucleotide_string(get(i)))
            .unwrap_or_default();

        let mismatches: u32 = col_mm.and_then(|i| get(i).parse().ok()).unwrap_or(0);

        // CFD and Hsu2013 are parsed from golden data but may be NA when scores
        // are disabled in production (calculate_cfd_score returns None).
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
            search_backend: "flashfry-1.15".to_string(),
            cfd_score,
            hsu_score,
        });
    }

    Ok(sites)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Strand;

    /// Parser accepts the documented FlashFry 1.15 column layout (native column
    /// names: `contig`, `start`, `stop`, `target`, `orientation`,
    /// `Doench2016CFDScore`, `Hsu2013`).
    /// This is a unit-parser test — NOT an E2E oracle parity test.
    #[test]
    fn parser_accepts_flashfry_1_15_fixture() {
        // Synthetic fixture matching documented FlashFry 1.15 output schema.
        // Coordinate convention: start=0-based, stop=0-based exclusive.
        let tsv = "\
contig\tstart\tstop\ttarget\torientation\tnumberOfMismatches\tDoench2016CFDScore\tHsu2013
chr\t99\t119\tGAGTCCGAGCAGAAGAAGAANGG\tFWD\t0\t1.0\t100.0
chr\t499\t519\tGAGTCCGAGCAGAAGAAGAANGG\tRVS\t2\t0.45\t72.3
";
        let sites = parse_flashfry(tsv).expect("parse should succeed");
        assert_eq!(sites.len(), 2);

        // Coordinate convention: start 99 (0-based) → 100 (1-based)
        // stop 119 (exclusive) → end = 119 (1-based inclusive last nt)
        assert_eq!(
            sites[0].start, 100,
            "0-based start should be converted to 1-based"
        );
        assert_eq!(
            sites[0].end, 119,
            "0-based exclusive stop → 1-based inclusive end"
        );
        assert_eq!(sites[0].strand, Strand::Plus, "FWD → Plus");
        // CFD and Hsu parsed from golden data
        assert_eq!(sites[0].cfd_score, Some(1.0));
        assert_eq!(sites[0].hsu_score, Some(100.0));

        assert_eq!(sites[1].strand, Strand::Minus, "RVS → Minus");
        assert_eq!(sites[1].mismatches, 2);
        assert_eq!(sites[1].search_backend, "flashfry-1.15");
    }

    /// Strand parser correctly maps FlashFry orientation tokens.
    #[test]
    fn strand_parser_accepts_fwd_rvs() {
        assert_eq!("FWD".parse::<Strand>(), Ok(Strand::Plus));
        assert_eq!("RVS".parse::<Strand>(), Ok(Strand::Minus));
        // Legacy tokens still work
        assert_eq!("+".parse::<Strand>(), Ok(Strand::Plus));
        assert_eq!("-".parse::<Strand>(), Ok(Strand::Minus));
    }

    /// NA / empty score fields are parsed as None (not a parse error).
    #[test]
    fn parser_handles_na_scores() {
        let tsv = "\
contig\tstart\tstop\ttarget\torientation\tnumberOfMismatches\tDoench2016CFDScore\tHsu2013
chr\t0\t23\tGAGTCCGAGCAGAAGAAGAANGG\tFWD\t0\tNA\tNA
";
        let sites = parse_flashfry(tsv).expect("parse should succeed");
        assert_eq!(sites[0].cfd_score, None);
        assert_eq!(sites[0].hsu_score, None);
    }
}
