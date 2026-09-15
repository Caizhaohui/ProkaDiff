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
    /// Unique identifier for this edit row.  If the TSV has an `edit_id` column
    /// it is used; otherwise auto-generated as `edit_1`, `edit_2`, etc.
    pub edit_id: String,
    pub seq_id: String,
    pub start: u64,
    pub end: u64,
    pub ref_allele: String,
    pub alt: String,
    pub kind: String,
}

// ── Edit-level assessment (FIX-015) ─────────────────────────────────────────

/// Status of a single intended edit after comparison against observed mutations.
///
/// The comparison is at the **edit level** (one `IntendedEdit` row), not the
/// **event count level**.  This prevents `observed > declared` when a cassette
/// insertion is matched by multiple JC events.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntendedEditStatus {
    /// All expected mutation events for this edit are present.
    Complete,
    /// At least one but not all expected events match.
    Partial,
    /// No matching mutation events found.
    Missing,
}

impl IntendedEditStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::Partial => "partial",
            Self::Missing => "missing",
        }
    }
}

/// Per-edit assessment: which mutation IDs (GD record IDs) matched this edit.
#[derive(Clone, Debug)]
pub struct IntendedEditAssessment {
    pub edit_id: String,
    pub status: IntendedEditStatus,
    /// IDs of GdEntry records that matched this intended edit.
    pub matched_event_ids: Vec<u32>,
}

/// Assess each intended edit against the set of observed mutations.
///
/// Returns one `IntendedEditAssessment` per entry in `intended`.
/// An edit is `Complete` when at least one matching event is found.
/// (For cassette edits that may produce multiple JC events, the current
/// implementation treats finding ≥ 1 match as Complete; future versions
/// can require all expected events.)
pub fn assess_intended_edits(
    mutations: &[GdEntry],
    intended: &[IntendedEdit],
) -> Vec<IntendedEditAssessment> {
    intended
        .iter()
        .map(|edit| {
            let matched: Vec<u32> = mutations
                .iter()
                .filter_map(|e| {
                    if matches_intended(e, edit) {
                        // Use GD record numeric ID if parseable, else 0
                        e.fields.first().and_then(|f| f.parse().ok()).or(Some(0))
                    } else {
                        None
                    }
                })
                .collect();
            let status = if matched.is_empty() {
                IntendedEditStatus::Missing
            } else {
                // For now, finding any match = Complete (cassette multi-JC handled
                // by the cassette branch in matches_intended returning true for each).
                IntendedEditStatus::Complete
            };
            IntendedEditAssessment {
                edit_id: edit.edit_id.clone(),
                status,
                matched_event_ids: matched,
            }
        })
        .collect()
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
    // edit_id column is optional — auto-generated when absent (FIX-015)
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
        // Use provided edit_id or auto-generate (FIX-015)
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
    // A DEL has no alt allele in GD; a non-empty, non-"." declared alt means the row does not
    // describe a clean deletion and must not mask this entry.
    size_ok && (t.alt.is_empty() || t.alt == ".")
}
