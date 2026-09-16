use prokadiff_gd::{GdEntry, GdKind};

#[derive(Debug, thiserror::Error)]
pub enum IntendedError {
    #[error("intended.tsv: {0}")]
    Parse(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntendedEdit {
    /// Unique identifier for this edit row. If the TSV has an `edit_id` column
    /// it is used; otherwise auto-generated as `edit_1`, `edit_2`, etc.
    pub edit_id: String,
    pub seq_id: String,
    pub start: u64,
    pub end: u64,
    pub ref_allele: String,
    pub alt: String,
    pub kind: String,
}

// ── Edit-level assessment ─────────────────────────────────────────

/// Status of a single intended edit after comprehensive evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntendedEditStatus {
    /// All expected mutation events/junctions for this edit are present and verified.
    Complete,
    /// At least one expected event/junction matches, but required counterpart or size is incomplete.
    Partial,
    /// No matching mutation events found at this locus.
    Missing,
    /// Target locus exhibits aberrant rearrangement, secondary junction, or unexpected insertion.
    UnexpectedStructure,
}

impl IntendedEditStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::Partial => "partial",
            Self::Missing => "missing",
            Self::UnexpectedStructure => "unexpected_structure",
        }
    }
}

/// Boundary assessment for structural edits (deletions, cassettes).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundaryAssessment {
    pub expected_pos: u64,
    pub observed_pos: Option<u64>,
    pub diff_bp: i64,
    pub passed: bool,
}

/// Per-edit comprehensive assessment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntendedEditAssessment {
    pub edit_id: String,
    pub kind: String,
    pub seq_id: String,
    pub expected_start: u64,
    pub expected_end: u64,
    pub status: IntendedEditStatus,
    /// IDs of GdEntry records that matched this intended edit (using real GdEntry.id).
    pub matched_event_ids: Vec<u32>,
    /// Optional boundary assessments for structural edits.
    pub left_boundary: Option<BoundaryAssessment>,
    pub right_boundary: Option<BoundaryAssessment>,
    /// Expected vs observed event size (e.g. deletion length).
    pub expected_size: Option<u64>,
    pub observed_size: Option<u64>,
    /// IDs of unexpected/aberrant secondary events at or near this locus.
    pub unexpected_event_ids: Vec<u32>,
    /// Human and machine readable notes.
    pub notes: Vec<String>,
}

pub fn parse_intended(text: &str) -> Result<Vec<IntendedEdit>, IntendedError> {
    let mut lines = text.lines().filter(|l| {
        let t = l.trim();
        !t.is_empty() && !t.starts_with('#')
    });
    let header = lines
        .next()
        .ok_or_else(|| IntendedError::Parse("missing header row".into()))?;
    let cols: Vec<String> = header
        .split('\t')
        .map(|s| s.trim().to_ascii_lowercase())
        .collect();
    let idx_or_alias = |primary: &str, alias: &str| -> Result<usize, IntendedError> {
        if let Some(pos) = cols.iter().position(|c| c == primary) {
            Ok(pos)
        } else if let Some(pos) = cols.iter().position(|c| c == alias) {
            eprintln!(
                "warning: deprecated column name '{alias}' in intended table; please use '{primary}'"
            );
            Ok(pos)
        } else {
            Err(IntendedError::Parse(format!(
                "missing column {primary} (or legacy alias {alias})"
            )))
        }
    };
    let idx = |name: &str| -> Result<usize, IntendedError> {
        cols.iter()
            .position(|c| c == name)
            .ok_or_else(|| IntendedError::Parse(format!("missing column {name}")))
    };
    let i_edit_id: Option<usize> = cols.iter().position(|c| c == "edit_id");
    let i_seq = idx("seq_id")?;
    let i_start = idx_or_alias("start", "position")?;
    let i_end = idx("end")?;
    let i_ref = idx("ref")?;
    let i_alt = idx("alt")?;
    let i_kind = idx_or_alias("kind", "gd_type")?;
    let mut out = Vec::new();
    for (n, line) in lines.enumerate() {
        let parts: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| parts.get(i).copied().unwrap_or("").trim();
        let start: u64 = get(i_start)
            .parse()
            .map_err(|_| IntendedError::Parse(format!("line {}: bad start", n + 2)))?;
        let end: u64 = get(i_end)
            .parse()
            .map_err(|_| IntendedError::Parse(format!("line {}: bad end", n + 2)))?;
        let edit_id = match i_edit_id {
            Some(idx) => {
                let s = get(idx);
                if s.is_empty() {
                    format!("edit_{}", n + 1)
                } else {
                    s.to_string()
                }
            }
            None => format!("edit_{}", n + 1),
        };
        out.push(IntendedEdit {
            edit_id,
            seq_id: get(i_seq).to_string(),
            start,
            end,
            ref_allele: get(i_ref).to_string(),
            alt: get(i_alt).to_string(),
            kind: get(i_kind).to_ascii_lowercase(),
        });
    }
    Ok(out)
}

pub fn parse_intended_path(
    path: impl AsRef<std::path::Path>,
) -> Result<Vec<IntendedEdit>, IntendedError> {
    let text = std::fs::read_to_string(path)?;
    parse_intended(&text)
}

/// Split mutations into (remaining, matched-intended).
pub fn mask_intended<'a>(
    mutations: &'a [GdEntry],
    intended: &[IntendedEdit],
) -> (Vec<&'a GdEntry>, Vec<&'a GdEntry>) {
    if intended.is_empty() {
        return (mutations.iter().collect(), Vec::new());
    }
    let mut remain = Vec::new();
    let mut observed = Vec::new();
    for e in mutations {
        if intended.iter().any(|t| matches_intended(e, t)) {
            observed.push(e);
        } else {
            remain.push(e);
        }
    }
    (remain, observed)
}

pub(crate) fn entry_intervals(e: &GdEntry) -> Vec<(String, u64, u64)> {
    match e.kind {
        GdKind::Ra | GdKind::Un | GdKind::Mc => Vec::new(),
        GdKind::Jc => {
            let mut v = Vec::new();
            if let (Some(s), Some(Ok(p))) =
                (e.fields.first(), e.fields.get(1).map(|x| x.parse::<u64>()))
            {
                v.push((s.clone(), p, p));
            }
            if let (Some(s), Some(Ok(p))) =
                (e.fields.get(3), e.fields.get(4).map(|x| x.parse::<u64>()))
            {
                v.push((s.clone(), p, p));
            }
            v
        }
        GdKind::Del | GdKind::Sub => {
            let Some(seq) = e.fields.first() else {
                return Vec::new();
            };
            let Some(Ok(start)) = e.fields.get(1).map(|x| x.parse::<u64>()) else {
                return Vec::new();
            };
            let size = e
                .fields
                .get(2)
                .and_then(|x| x.parse::<u64>().ok())
                .unwrap_or(1)
                .max(1);
            vec![(seq.clone(), start, start.saturating_add(size - 1))]
        }
        _ => match (e.seq_id(), e.position()) {
            (Some(s), Some(p)) => vec![(s.to_string(), p, p)],
            _ => Vec::new(),
        },
    }
}

fn allele(e: &GdEntry) -> Option<&str> {
    match e.kind {
        GdKind::Snp | GdKind::Ins => e.fields.get(2).map(String::as_str),
        GdKind::Sub => e.fields.get(3).map(String::as_str),
        _ => None,
    }
}

fn alt_ok(e: &GdEntry, t: &IntendedEdit) -> bool {
    if t.alt.is_empty() || t.alt == "." {
        return true;
    }
    match allele(e) {
        Some(a) => a.eq_ignore_ascii_case(&t.alt),
        None => true,
    }
}

fn overlaps_intended(e: &GdEntry, t: &IntendedEdit) -> bool {
    entry_intervals(e)
        .iter()
        .any(|(sid, a, b)| sid == &t.seq_id && *a <= t.end && t.start <= *b)
}

fn is_near_cassette(e: &GdEntry, t: &IntendedEdit) -> bool {
    let window = 50;
    let w_start = t.start.saturating_sub(window);
    let w_end = t.end.saturating_add(window);
    entry_intervals(e)
        .iter()
        .any(|(sid, a, b)| sid == &t.seq_id && *a <= w_end && w_start <= *b)
}

fn matches_intended(e: &GdEntry, t: &IntendedEdit) -> bool {
    if !overlaps_intended(e, t) {
        return false;
    }
    match t.kind.as_str() {
        "snp" => e.kind == GdKind::Snp && alt_ok(e, t),
        "sub" => e.kind == GdKind::Sub && alt_ok(e, t),
        "ins" => e.kind == GdKind::Ins && alt_ok(e, t),
        "del" | "indel" => del_matches_intended(e, t),
        "cassette" => {
            matches!(e.kind, GdKind::Jc | GdKind::Mob | GdKind::Ins | GdKind::Con)
        }
        _ => false,
    }
}

fn del_size(e: &GdEntry) -> Option<u64> {
    e.fields.get(2).and_then(|s| s.parse().ok())
}

fn del_matches_intended(e: &GdEntry, t: &IntendedEdit) -> bool {
    if e.kind != GdKind::Del {
        return false;
    }
    let span = t.end.saturating_sub(t.start) + 1;
    let size_ok = match del_size(e) {
        Some(s) => s == span,
        None => false,
    };
    size_ok && (t.alt.is_empty() || t.alt == ".")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CassetteJunctionSide {
    Left,
    Right,
    Ambiguous,
}

fn classify_cassette_junction(
    e: &GdEntry,
    edit: &IntendedEdit,
    window: u64,
) -> Option<(CassetteJunctionSide, u64, i64)> {
    let mut hits = Vec::new();
    if e.kind == GdKind::Jc {
        if let (Some(s1), Some(Ok(p1)), Some(st1)) = (
            e.fields.first(),
            e.fields.get(1).map(|x| x.parse::<u64>()),
            e.fields.get(2),
        ) {
            if s1 == &edit.seq_id {
                hits.push((p1, st1.as_str()));
            }
        }
        if let (Some(s2), Some(Ok(p2)), Some(st2)) = (
            e.fields.get(3),
            e.fields.get(4).map(|x| x.parse::<u64>()),
            e.fields.get(5),
        ) {
            if s2 == &edit.seq_id {
                hits.push((p2, st2.as_str()));
            }
        }
    } else if let (Some(s), Some(p)) = (e.seq_id(), e.position()) {
        if s == edit.seq_id {
            hits.push((p, "+"));
        }
    }

    let mut best: Option<(CassetteJunctionSide, u64, i64, u64)> = None;
    for (pos, strand) in hits {
        let d_start = (pos as i64 - edit.start as i64).unsigned_abs();
        let d_end = (pos as i64 - edit.end as i64).unsigned_abs();

        if d_start > window && d_end > window {
            continue;
        }

        let (side, diff, dist) = if edit.start == edit.end
            || (edit.end.saturating_sub(edit.start) <= 5 && d_start <= window && d_end <= window)
        {
            // Insertion at single locus: use strand to differentiate flanks
            if strand == "+" || strand == "0" || strand == "F" || strand == "forward" {
                (
                    CassetteJunctionSide::Left,
                    pos as i64 - edit.start as i64,
                    d_start,
                )
            } else if strand == "-" || strand == "1" || strand == "R" || strand == "reverse" {
                (
                    CassetteJunctionSide::Right,
                    pos as i64 - edit.end as i64,
                    d_end,
                )
            } else {
                (
                    CassetteJunctionSide::Ambiguous,
                    pos as i64 - edit.start as i64,
                    d_start,
                )
            }
        } else if d_start < d_end && d_start <= window {
            (
                CassetteJunctionSide::Left,
                pos as i64 - edit.start as i64,
                d_start,
            )
        } else if d_end < d_start && d_end <= window {
            (
                CassetteJunctionSide::Right,
                pos as i64 - edit.end as i64,
                d_end,
            )
        } else if d_start <= window {
            (
                CassetteJunctionSide::Ambiguous,
                pos as i64 - edit.start as i64,
                d_start,
            )
        } else {
            continue;
        };

        match best {
            None => best = Some((side, pos, diff, dist)),
            Some((_, _, _, prev_dist)) if dist < prev_dist => {
                best = Some((side, pos, diff, dist));
            }
            _ => {}
        }
    }

    best.map(|(side, pos, diff, _)| (side, pos, diff))
}

/// Assess a single intended edit against the set of observed mutations.
fn assess_single_edit(edit: &IntendedEdit, mutations: &[GdEntry]) -> IntendedEditAssessment {
    let kind = edit.kind.to_ascii_lowercase();
    let expected_span = edit.end.saturating_sub(edit.start) + 1;

    let mut matched_event_ids = Vec::new();
    let mut unexpected_event_ids = Vec::new();
    let mut notes = Vec::new();
    let mut left_boundary = None;
    let mut right_boundary = None;
    let mut expected_size = None;
    let mut observed_size = None;

    match kind.as_str() {
        "snp" | "sub" | "ins" => {
            let target_kind = match kind.as_str() {
                "snp" => GdKind::Snp,
                "sub" => GdKind::Sub,
                "ins" => GdKind::Ins,
                _ => unreachable!(),
            };
            let mut exact_matches = Vec::new();
            let mut allele_mismatches = Vec::new();
            let mut other_events = Vec::new();

            for e in mutations {
                if overlaps_intended(e, edit) {
                    if e.kind == target_kind {
                        if alt_ok(e, edit) {
                            exact_matches.push(e);
                        } else {
                            allele_mismatches.push(e);
                        }
                    } else if matches!(e.kind, GdKind::Jc | GdKind::Del | GdKind::Mob) {
                        other_events.push(e);
                    }
                }
            }

            let status = if !exact_matches.is_empty() {
                for e in &exact_matches {
                    matched_event_ids.push(e.id);
                }
                if !other_events.is_empty() {
                    for e in &other_events {
                        unexpected_event_ids.push(e.id);
                    }
                    notes.push("additional unexpected structural event at target locus".into());
                    IntendedEditStatus::UnexpectedStructure
                } else {
                    IntendedEditStatus::Complete
                }
            } else if !allele_mismatches.is_empty() {
                for e in &allele_mismatches {
                    matched_event_ids.push(e.id);
                    if let Some(obs) = allele(e) {
                        notes.push(format!(
                            "mismatched allele: expected {}, observed {}",
                            edit.alt, obs
                        ));
                    }
                }
                IntendedEditStatus::Partial
            } else if !other_events.is_empty() {
                for e in &other_events {
                    unexpected_event_ids.push(e.id);
                }
                notes.push("target site disrupted by unexpected structural variant".into());
                IntendedEditStatus::UnexpectedStructure
            } else {
                IntendedEditStatus::Missing
            };

            IntendedEditAssessment {
                edit_id: edit.edit_id.clone(),
                kind: edit.kind.clone(),
                seq_id: edit.seq_id.clone(),
                expected_start: edit.start,
                expected_end: edit.end,
                status,
                matched_event_ids,
                left_boundary,
                right_boundary,
                expected_size,
                observed_size,
                unexpected_event_ids,
                notes,
            }
        }
        "del" | "indel" => {
            expected_size = Some(expected_span);
            let mut del_matches = Vec::new();
            let mut jc_matches = Vec::new();
            let mut other_events = Vec::new();

            for e in mutations {
                if overlaps_intended(e, edit) {
                    if e.kind == GdKind::Del {
                        del_matches.push(e);
                    } else if e.kind == GdKind::Jc {
                        jc_matches.push(e);
                    } else {
                        other_events.push(e);
                    }
                }
            }

            let status = if let Some(e) = del_matches.first() {
                matched_event_ids.push(e.id);
                let d_size = del_size(e).unwrap_or(1);
                observed_size = Some(d_size);
                let obs_start = e.position().unwrap_or(edit.start);
                let obs_end = obs_start.saturating_add(d_size.saturating_sub(1));

                let left_diff = (obs_start as i64) - (edit.start as i64);
                let right_diff = (obs_end as i64) - (edit.end as i64);
                let left_pass = left_diff.abs() <= 2;
                let right_pass = right_diff.abs() <= 2;

                left_boundary = Some(BoundaryAssessment {
                    expected_pos: edit.start,
                    observed_pos: Some(obs_start),
                    diff_bp: left_diff,
                    passed: left_pass,
                });
                right_boundary = Some(BoundaryAssessment {
                    expected_pos: edit.end,
                    observed_pos: Some(obs_end),
                    diff_bp: right_diff,
                    passed: right_pass,
                });

                if left_pass && right_pass && (edit.alt.is_empty() || edit.alt == ".") {
                    IntendedEditStatus::Complete
                } else {
                    notes.push(format!(
                        "boundary or size discrepancy: left diff {} bp, right diff {} bp",
                        left_diff, right_diff
                    ));
                    IntendedEditStatus::Partial
                }
            } else if !jc_matches.is_empty() {
                for jc in &jc_matches {
                    matched_event_ids.push(jc.id);
                }
                notes.push("deletion supported by junction evidence".into());
                IntendedEditStatus::Complete
            } else if !other_events.is_empty() {
                for e in &other_events {
                    unexpected_event_ids.push(e.id);
                }
                notes.push("unexpected variant at deletion locus".into());
                IntendedEditStatus::UnexpectedStructure
            } else {
                IntendedEditStatus::Missing
            };

            IntendedEditAssessment {
                edit_id: edit.edit_id.clone(),
                kind: edit.kind.clone(),
                seq_id: edit.seq_id.clone(),
                expected_start: edit.start,
                expected_end: edit.end,
                status,
                matched_event_ids,
                left_boundary,
                right_boundary,
                expected_size,
                observed_size,
                unexpected_event_ids,
                notes,
            }
        }
        "cassette" => {
            let mut jc_matches = Vec::new();
            let mut other_events = Vec::new();

            for e in mutations {
                if is_near_cassette(e, edit) {
                    if matches!(e.kind, GdKind::Jc | GdKind::Mob) {
                        jc_matches.push(e);
                    } else {
                        other_events.push(e);
                    }
                }
            }

            let mut left_candidates = Vec::new();
            let mut right_candidates = Vec::new();
            let mut ambiguous_candidates = Vec::new();

            for e in &jc_matches {
                if let Some((side, pos, diff)) = classify_cassette_junction(e, edit, 50) {
                    match side {
                        CassetteJunctionSide::Left => left_candidates.push((e, pos, diff)),
                        CassetteJunctionSide::Right => right_candidates.push((e, pos, diff)),
                        CassetteJunctionSide::Ambiguous => {
                            ambiguous_candidates.push((e, pos, diff))
                        }
                    }
                }
            }

            // Differentiate ambiguous candidates when one side is missing
            if left_candidates.is_empty()
                && !ambiguous_candidates.is_empty()
                && right_candidates.len() == 1
            {
                let cand = ambiguous_candidates.remove(0);
                left_candidates.push(cand);
            } else if right_candidates.is_empty()
                && !ambiguous_candidates.is_empty()
                && left_candidates.len() == 1
            {
                let cand = ambiguous_candidates.remove(0);
                right_candidates.push(cand);
            } else if left_candidates.is_empty()
                && right_candidates.is_empty()
                && ambiguous_candidates.len() == 2
            {
                let cand1 = ambiguous_candidates.remove(0);
                let cand2 = ambiguous_candidates.remove(0);
                left_candidates.push(cand1);
                right_candidates.push(cand2);
            }

            let status = if left_candidates.len() == 1
                && right_candidates.len() == 1
                && jc_matches.len() == 2
                && other_events.is_empty()
            {
                let (left_e, left_pos, left_diff) = left_candidates[0];
                let (right_e, right_pos, right_diff) = right_candidates[0];

                matched_event_ids.push(left_e.id);
                matched_event_ids.push(right_e.id);

                let left_pass = left_diff.abs() <= 5;
                let right_pass = right_diff.abs() <= 5;

                left_boundary = Some(BoundaryAssessment {
                    expected_pos: edit.start,
                    observed_pos: Some(left_pos),
                    diff_bp: left_diff,
                    passed: left_pass,
                });
                right_boundary = Some(BoundaryAssessment {
                    expected_pos: edit.end,
                    observed_pos: Some(right_pos),
                    diff_bp: right_diff,
                    passed: right_pass,
                });

                if left_pass && right_pass {
                    notes.push("both cassette junctions confirmed".into());
                    IntendedEditStatus::Complete
                } else {
                    notes.push(format!(
                        "cassette junction boundary discrepancy: left diff {} bp, right diff {} bp",
                        left_diff, right_diff
                    ));
                    IntendedEditStatus::Partial
                }
            } else if left_candidates.len() == 1
                && right_candidates.is_empty()
                && jc_matches.len() == 1
                && other_events.is_empty()
            {
                let (left_e, left_pos, left_diff) = left_candidates[0];
                matched_event_ids.push(left_e.id);
                left_boundary = Some(BoundaryAssessment {
                    expected_pos: edit.start,
                    observed_pos: Some(left_pos),
                    diff_bp: left_diff,
                    passed: left_diff.abs() <= 5,
                });
                notes.push(
                    "single junction detected (left flank candidate, partial integration)".into(),
                );
                IntendedEditStatus::Partial
            } else if right_candidates.len() == 1
                && left_candidates.is_empty()
                && jc_matches.len() == 1
                && other_events.is_empty()
            {
                let (right_e, right_pos, right_diff) = right_candidates[0];
                matched_event_ids.push(right_e.id);
                right_boundary = Some(BoundaryAssessment {
                    expected_pos: edit.end,
                    observed_pos: Some(right_pos),
                    diff_bp: right_diff,
                    passed: right_diff.abs() <= 5,
                });
                notes.push(
                    "single junction detected (right flank candidate, partial integration)".into(),
                );
                IntendedEditStatus::Partial
            } else if left_candidates.len() >= 2 && right_candidates.is_empty() {
                for (e, _, _) in &left_candidates {
                    unexpected_event_ids.push(e.id);
                }
                notes.push("aberrant multiple left junctions detected at cassette locus".into());
                IntendedEditStatus::UnexpectedStructure
            } else if right_candidates.len() >= 2 && left_candidates.is_empty() {
                for (e, _, _) in &right_candidates {
                    unexpected_event_ids.push(e.id);
                }
                notes.push("aberrant multiple right junctions detected at cassette locus".into());
                IntendedEditStatus::UnexpectedStructure
            } else if left_candidates.len() == 1 && right_candidates.len() == 1 {
                // Correct pair + extra events
                let (left_e, left_pos, left_diff) = left_candidates[0];
                let (right_e, right_pos, right_diff) = right_candidates[0];
                matched_event_ids.push(left_e.id);
                matched_event_ids.push(right_e.id);
                left_boundary = Some(BoundaryAssessment {
                    expected_pos: edit.start,
                    observed_pos: Some(left_pos),
                    diff_bp: left_diff,
                    passed: left_diff.abs() <= 5,
                });
                right_boundary = Some(BoundaryAssessment {
                    expected_pos: edit.end,
                    observed_pos: Some(right_pos),
                    diff_bp: right_diff,
                    passed: right_diff.abs() <= 5,
                });
                for e in &jc_matches {
                    if e.id != left_e.id && e.id != right_e.id {
                        unexpected_event_ids.push(e.id);
                    }
                }
                for e in &other_events {
                    unexpected_event_ids.push(e.id);
                }
                notes.push("correct junction pair detected with additional aberrant events at cassette locus".into());
                IntendedEditStatus::UnexpectedStructure
            } else if !jc_matches.is_empty() || !other_events.is_empty() {
                for e in &jc_matches {
                    unexpected_event_ids.push(e.id);
                }
                for e in &other_events {
                    unexpected_event_ids.push(e.id);
                }
                let msg = if jc_matches.len() > 1 {
                    "aberrant multiple junctions detected at cassette locus"
                } else {
                    "aberrant or unrecognized junction arrangement at cassette locus"
                };
                notes.push(msg.into());
                IntendedEditStatus::UnexpectedStructure
            } else {
                IntendedEditStatus::Missing
            };

            IntendedEditAssessment {
                edit_id: edit.edit_id.clone(),
                kind: edit.kind.clone(),
                seq_id: edit.seq_id.clone(),
                expected_start: edit.start,
                expected_end: edit.end,
                status,
                matched_event_ids,
                left_boundary,
                right_boundary,
                expected_size,
                observed_size,
                unexpected_event_ids,
                notes,
            }
        }
        _ => {
            let matched: Vec<u32> = mutations
                .iter()
                .filter(|e| matches_intended(e, edit))
                .map(|e| e.id)
                .collect();
            let status = if matched.is_empty() {
                IntendedEditStatus::Missing
            } else {
                IntendedEditStatus::Complete
            };
            IntendedEditAssessment {
                edit_id: edit.edit_id.clone(),
                kind: edit.kind.clone(),
                seq_id: edit.seq_id.clone(),
                expected_start: edit.start,
                expected_end: edit.end,
                status,
                matched_event_ids: matched,
                left_boundary: None,
                right_boundary: None,
                expected_size: None,
                observed_size: None,
                unexpected_event_ids: Vec::new(),
                notes: Vec::new(),
            }
        }
    }
}

/// Assess each intended edit against the set of observed mutations.
///
/// Returns one `IntendedEditAssessment` per entry in `intended`.
pub fn assess_intended_edits(
    mutations: &[GdEntry],
    intended: &[IntendedEdit],
) -> Vec<IntendedEditAssessment> {
    intended
        .iter()
        .map(|edit| assess_single_edit(edit, mutations))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_matched_event_ids_uses_real_gd_id() {
        let entry = GdEntry::snp(42, "NC_000913.3", 100, "T");
        let edit = IntendedEdit {
            edit_id: "edit_1".into(),
            seq_id: "NC_000913.3".into(),
            start: 100,
            end: 100,
            ref_allele: "A".into(),
            alt: "T".into(),
            kind: "snp".into(),
        };
        let assessments = assess_intended_edits(&[entry], &[edit]);
        assert_eq!(assessments.len(), 1);
        assert_eq!(assessments[0].status, IntendedEditStatus::Complete);
        assert_eq!(assessments[0].matched_event_ids, vec![42]);
    }

    #[test]
    fn test_snp_allele_mismatch_yields_partial() {
        let entry = GdEntry::snp(10, "chr", 100, "C");
        let edit = IntendedEdit {
            edit_id: "edit_snp".into(),
            seq_id: "chr".into(),
            start: 100,
            end: 100,
            ref_allele: "A".into(),
            alt: "T".into(),
            kind: "snp".into(),
        };
        let assessments = assess_intended_edits(&[entry], &[edit]);
        assert_eq!(assessments[0].status, IntendedEditStatus::Partial);
        assert_eq!(assessments[0].matched_event_ids, vec![10]);
        assert!(assessments[0].notes[0].contains("mismatched allele"));
    }

    #[test]
    fn test_large_del_boundary_assessment_complete() {
        let entry = GdEntry::del(15, "chr", 1000, 500); // 1000..=1499
        let edit = IntendedEdit {
            edit_id: "edit_del".into(),
            seq_id: "chr".into(),
            start: 1000,
            end: 1499,
            ref_allele: ".".into(),
            alt: ".".into(),
            kind: "del".into(),
        };
        let assessments = assess_intended_edits(&[entry], &[edit]);
        assert_eq!(assessments[0].status, IntendedEditStatus::Complete);
        assert_eq!(assessments[0].matched_event_ids, vec![15]);
        assert_eq!(assessments[0].expected_size, Some(500));
        assert_eq!(assessments[0].observed_size, Some(500));
        assert!(assessments[0].left_boundary.as_ref().unwrap().passed);
        assert!(assessments[0].right_boundary.as_ref().unwrap().passed);
    }

    #[test]
    fn test_large_del_boundary_mismatch_yields_partial() {
        let entry = GdEntry::del(16, "chr", 1000, 300); // 1000..=1299 vs declared 1000..=1499
        let edit = IntendedEdit {
            edit_id: "edit_del_partial".into(),
            seq_id: "chr".into(),
            start: 1000,
            end: 1499,
            ref_allele: ".".into(),
            alt: ".".into(),
            kind: "del".into(),
        };
        let assessments = assess_intended_edits(&[entry], &[edit]);
        assert_eq!(assessments[0].status, IntendedEditStatus::Partial);
        assert_eq!(assessments[0].matched_event_ids, vec![16]);
        assert_eq!(assessments[0].observed_size, Some(300));
        assert!(!assessments[0].right_boundary.as_ref().unwrap().passed);
    }

    #[test]
    fn test_cassette_single_jc_yields_partial() {
        let jc1 = GdEntry::jc(21, "chr", 5000, "+", "plasmid", 100, "-", 0);
        let edit = IntendedEdit {
            edit_id: "cassette_1".into(),
            seq_id: "chr".into(),
            start: 5000,
            end: 5000,
            ref_allele: ".".into(),
            alt: ".".into(),
            kind: "cassette".into(),
        };
        let assessments = assess_intended_edits(&[jc1], &[edit]);
        assert_eq!(assessments[0].status, IntendedEditStatus::Partial);
        assert_eq!(assessments[0].matched_event_ids, vec![21]);
        assert!(assessments[0].notes[0].contains("single junction detected"));
    }

    #[test]
    fn test_cassette_two_jc_yields_complete() {
        let jc1 = GdEntry::jc(31, "chr", 5000, "+", "plasmid", 100, "-", 0);
        let jc2 = GdEntry::jc(32, "chr", 5000, "-", "plasmid", 2500, "+", 0);
        let edit = IntendedEdit {
            edit_id: "cassette_1".into(),
            seq_id: "chr".into(),
            start: 5000,
            end: 5000,
            ref_allele: ".".into(),
            alt: ".".into(),
            kind: "cassette".into(),
        };
        let assessments = assess_intended_edits(&[jc1, jc2], &[edit]);
        assert_eq!(assessments[0].status, IntendedEditStatus::Complete);
        assert_eq!(assessments[0].matched_event_ids, vec![31, 32]);
        assert!(assessments[0].left_boundary.as_ref().unwrap().passed);
        assert!(assessments[0].right_boundary.as_ref().unwrap().passed);
    }

    #[test]
    fn test_cassette_three_jc_yields_unexpected_structure() {
        let jc1 = GdEntry::jc(41, "chr", 5000, "+", "plasmid", 100, "-", 0);
        let jc2 = GdEntry::jc(42, "chr", 5000, "-", "plasmid", 2500, "+", 0);
        let jc3 = GdEntry::jc(43, "chr", 5010, "+", "plasmid", 500, "+", 0);
        let edit = IntendedEdit {
            edit_id: "cassette_aberrant".into(),
            seq_id: "chr".into(),
            start: 5000,
            end: 5000,
            ref_allele: ".".into(),
            alt: ".".into(),
            kind: "cassette".into(),
        };
        let assessments = assess_intended_edits(&[jc1, jc2, jc3], &[edit]);
        assert_eq!(
            assessments[0].status,
            IntendedEditStatus::UnexpectedStructure
        );
        assert!(assessments[0].notes[0].contains("aberrant multiple junctions"));
    }

    #[test]
    fn test_cassette_two_left_jc_yields_unexpected_structure() {
        // Both JCs have positive strand at the same position -> 2 left-like JCs
        let jc1 = GdEntry::jc(51, "chr", 5000, "+", "plasmid", 100, "-", 0);
        let jc2 = GdEntry::jc(52, "chr", 5000, "+", "plasmid", 2500, "+", 0);
        let edit = IntendedEdit {
            edit_id: "cassette_two_left".into(),
            seq_id: "chr".into(),
            start: 5000,
            end: 5000,
            ref_allele: ".".into(),
            alt: ".".into(),
            kind: "cassette".into(),
        };
        let assessments = assess_intended_edits(&[jc1, jc2], &[edit]);
        assert_eq!(
            assessments[0].status,
            IntendedEditStatus::UnexpectedStructure
        );
        assert!(assessments[0].notes[0].contains("aberrant multiple left junctions"));
    }

    #[test]
    fn test_cassette_replacement_uses_actual_coordinates() {
        // Replacement spanning 1000..=2000
        // Left junction at 1002 (diff +2), Right junction at 1999 (diff -1)
        let jc_left = GdEntry::jc(61, "chr", 1002, "+", "donor", 1, "+", 0);
        let jc_right = GdEntry::jc(62, "chr", 1999, "-", "donor", 500, "-", 0);
        let edit = IntendedEdit {
            edit_id: "cassette_rep".into(),
            seq_id: "chr".into(),
            start: 1000,
            end: 2000,
            ref_allele: ".".into(),
            alt: ".".into(),
            kind: "cassette".into(),
        };
        let assessments = assess_intended_edits(&[jc_left, jc_right], &[edit]);
        assert_eq!(assessments[0].status, IntendedEditStatus::Complete);
        let left_b = assessments[0].left_boundary.as_ref().unwrap();
        assert_eq!(left_b.expected_pos, 1000);
        assert_eq!(left_b.observed_pos, Some(1002));
        assert_eq!(left_b.diff_bp, 2);
        assert!(left_b.passed);

        let right_b = assessments[0].right_boundary.as_ref().unwrap();
        assert_eq!(right_b.expected_pos, 2000);
        assert_eq!(right_b.observed_pos, Some(1999));
        assert_eq!(right_b.diff_bp, -1);
        assert!(right_b.passed);
    }

    #[test]
    fn test_missing_intended_edit() {
        let edit = IntendedEdit {
            edit_id: "edit_missing".into(),
            seq_id: "chr".into(),
            start: 1000,
            end: 1000,
            ref_allele: "C".into(),
            alt: "G".into(),
            kind: "snp".into(),
        };
        let assessments = assess_intended_edits(&[], &[edit]);
        assert_eq!(assessments[0].status, IntendedEditStatus::Missing);
        assert!(assessments[0].matched_event_ids.is_empty());
    }
}
