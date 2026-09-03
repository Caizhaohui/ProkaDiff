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
    pub seq_id: String,
    pub start: u64,
    pub end: u64,
    pub ref_allele: String,
    pub alt: String,
    pub kind: String,
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
    let idx = |name: &str| -> Result<usize, IntendedError> {
        cols.iter()
            .position(|c| c == name)
            .ok_or_else(|| IntendedError::Parse(format!("missing column {name}")))
    };
    let i_seq = idx("seq_id")?;
    let i_start = idx("start")?;
    let i_end = idx("end")?;
    let i_ref = idx("ref")?;
    let i_alt = idx("alt")?;
    let i_kind = idx("kind")?;
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
        out.push(IntendedEdit {
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
        GdKind::Del => {
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
