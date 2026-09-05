//! Genome Diff subset I/O and `gdtools SUBTRACT`-style set difference.
//!
//! Matching key: mutation type + coordinates + allele (see `docs/schema.md`).

#![deny(unsafe_code)]

use std::collections::{BTreeMap, HashSet};
use std::fmt;
use std::fs;
use std::path::Path;
use std::str::FromStr;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum GdError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error on line {line}: {msg}")]
    Parse { line: usize, msg: String },
}

pub type Result<T> = std::result::Result<T, GdError>;

/// Three-letter mutation and two-letter evidence types in the first-period subset.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum GdKind {
    Snp,
    Ins,
    Del,
    Mob,
    Amp,
    Con,
    Jc,
    Un,
    Ra,
    Mc,
}

impl GdKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Snp => "SNP",
            Self::Ins => "INS",
            Self::Del => "DEL",
            Self::Mob => "MOB",
            Self::Amp => "AMP",
            Self::Con => "CON",
            Self::Jc => "JC",
            Self::Un => "UN",
            Self::Ra => "RA",
            Self::Mc => "MC",
        }
    }

    fn field_count(self) -> usize {
        match self {
            Self::Snp | Self::Ins => 3,
            Self::Del | Self::Un => 3,
            Self::Mob => 5,
            Self::Amp => 4,
            Self::Con => 4,
            Self::Jc => 7,
            Self::Ra => 5,
            Self::Mc => 5,
        }
    }
}

impl FromStr for GdKind {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        Ok(match s {
            "SNP" => Self::Snp,
            "INS" => Self::Ins,
            "DEL" => Self::Del,
            "MOB" => Self::Mob,
            "AMP" => Self::Amp,
            "CON" => Self::Con,
            "JC" => Self::Jc,
            "UN" => Self::Un,
            "RA" => Self::Ra,
            "MC" => Self::Mc,
            other => return Err(format!("unsupported GD type {other}")),
        })
    }
}

/// One Genome Diff record (mutation or evidence).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GdEntry {
    pub kind: GdKind,
    pub id: u32,
    pub parent_ids: Vec<u32>,
    /// Required positional fields in published GD column order.
    pub fields: Vec<String>,
    pub attrs: BTreeMap<String, String>,
}

impl GdEntry {
    pub fn snp(
        id: u32,
        seq_id: impl Into<String>,
        position: u64,
        new_seq: impl Into<String>,
    ) -> Self {
        Self {
            kind: GdKind::Snp,
            id,
            parent_ids: Vec::new(),
            fields: vec![seq_id.into(), position.to_string(), new_seq.into()],
            attrs: BTreeMap::new(),
        }
    }

    pub fn ins(
        id: u32,
        seq_id: impl Into<String>,
        position: u64,
        new_seq: impl Into<String>,
    ) -> Self {
        Self {
            kind: GdKind::Ins,
            id,
            parent_ids: Vec::new(),
            fields: vec![seq_id.into(), position.to_string(), new_seq.into()],
            attrs: BTreeMap::new(),
        }
    }

    pub fn del(id: u32, seq_id: impl Into<String>, position: u64, size: u64) -> Self {
        Self {
            kind: GdKind::Del,
            id,
            parent_ids: Vec::new(),
            fields: vec![seq_id.into(), position.to_string(), size.to_string()],
            attrs: BTreeMap::new(),
        }
    }

    pub fn mob(
        id: u32,
        seq_id: impl Into<String>,
        position: u64,
        repeat_name: impl Into<String>,
        strand: impl Into<String>,
        duplication_size: i64,
    ) -> Self {
        Self {
            kind: GdKind::Mob,
            id,
            parent_ids: Vec::new(),
            fields: vec![
                seq_id.into(),
                position.to_string(),
                repeat_name.into(),
                strand.into(),
                duplication_size.to_string(),
            ],
            attrs: BTreeMap::new(),
        }
    }

    pub fn mc(
        id: u32,
        seq_id: impl Into<String>,
        start: u64,
        end: u64,
        start_range: u64,
        end_range: u64,
    ) -> Self {
        Self {
            kind: GdKind::Mc,
            id,
            parent_ids: Vec::new(),
            fields: vec![
                seq_id.into(),
                start.to_string(),
                end.to_string(),
                start_range.to_string(),
                end_range.to_string(),
            ],
            attrs: BTreeMap::new(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn jc(
        id: u32,
        side_1_seq_id: impl Into<String>,
        side_1_position: u64,
        side_1_strand: impl Into<String>,
        side_2_seq_id: impl Into<String>,
        side_2_position: u64,
        side_2_strand: impl Into<String>,
        overlap: i64,
    ) -> Self {
        Self {
            kind: GdKind::Jc,
            id,
            parent_ids: Vec::new(),
            fields: vec![
                side_1_seq_id.into(),
                side_1_position.to_string(),
                side_1_strand.into(),
                side_2_seq_id.into(),
                side_2_position.to_string(),
                side_2_strand.into(),
                overlap.to_string(),
            ],
            attrs: BTreeMap::new(),
        }
    }

    pub fn seq_id(&self) -> Option<&str> {
        match self.kind {
            GdKind::Jc => self.fields.first().map(String::as_str),
            _ => self.fields.first().map(String::as_str),
        }
    }

    pub fn position(&self) -> Option<u64> {
        match self.kind {
            GdKind::Mc | GdKind::Un => self.fields.get(1).and_then(|s| s.parse().ok()),
            _ => self.fields.get(1).and_then(|s| s.parse().ok()),
        }
    }

    pub fn del_size(&self) -> Option<u64> {
        if self.kind == GdKind::Del {
            self.fields.get(2).and_then(|s| s.parse().ok())
        } else {
            None
        }
    }

    /// Set-difference key: type + coordinates + allele.
    pub fn subtract_key(&self) -> String {
        format!("{}|{}", self.kind.as_str(), self.fields.join("|"))
    }
}

/// A Genome Diff document.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GenomeDiff {
    pub metadata: Vec<(String, String)>,
    pub entries: Vec<GdEntry>,
}

impl GenomeDiff {
    pub fn new() -> Self {
        Self {
            metadata: vec![("GENOME_DIFF".into(), "1.0".into())],
            entries: Vec::new(),
        }
    }

    pub fn parse(text: &str) -> Result<Self> {
        let mut doc = Self {
            metadata: Vec::new(),
            entries: Vec::new(),
        };
        for (i, raw) in text.lines().enumerate() {
            let line_no = i + 1;
            let line = raw.trim_end();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("#=") {
                let rest = rest.trim();
                let (k, v) = rest
                    .split_once(char::is_whitespace)
                    .map(|(k, v)| (k.to_string(), v.trim().to_string()))
                    .unwrap_or_else(|| (rest.to_string(), String::new()));
                doc.metadata.push((k, v));
                continue;
            }
            if line.starts_with('#') {
                continue;
            }
            doc.entries.push(parse_entry(line, line_no)?);
        }
        if !doc.metadata.iter().any(|(k, _)| k == "GENOME_DIFF") {
            doc.metadata.insert(0, ("GENOME_DIFF".into(), "1.0".into()));
        }
        Ok(doc)
    }

    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let text = fs::read_to_string(path)?;
        Self::parse(&text)
    }

    pub fn to_gd_string(&self) -> String {
        let mut out = String::new();
        for (k, v) in &self.metadata {
            if v.is_empty() {
                out.push_str(&format!("#={k}\n"));
            } else {
                out.push_str(&format!("#={k}\t{v}\n"));
            }
        }
        for e in &self.entries {
            out.push_str(&e.to_line());
            out.push('\n');
        }
        out
    }

    pub fn write_path(&self, path: impl AsRef<Path>) -> Result<()> {
        fs::write(path, self.to_gd_string())?;
        Ok(())
    }

    /// `gdtools SUBTRACT self other`: records in `self` whose subtract key is not in `other`.
    /// UN evidence is kept even when a matching UN exists in `other` (parity: UN is not an
    /// oracle failure). Mutation and other evidence types are subtracted.
    pub fn subtract(&self, other: &GenomeDiff) -> GenomeDiff {
        let remove: HashSet<String> = other
            .entries
            .iter()
            .filter(|e| e.kind != GdKind::Un)
            .map(GdEntry::subtract_key)
            .collect();
        GenomeDiff {
            metadata: self.metadata.clone(),
            entries: self
                .entries
                .iter()
                .filter(|e| e.kind == GdKind::Un || !remove.contains(&e.subtract_key()))
                .cloned()
                .collect(),
        }
    }
}

impl GdEntry {
    fn to_line(&self) -> String {
        let parents = if self.parent_ids.is_empty() {
            ".".to_string()
        } else {
            self.parent_ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };
        let mut cols = vec![self.kind.as_str().to_string(), self.id.to_string(), parents];
        cols.extend(self.fields.iter().cloned());
        for (k, v) in &self.attrs {
            cols.push(format!("{k}={v}"));
        }
        cols.join("\t")
    }
}

fn parse_entry(line: &str, line_no: usize) -> Result<GdEntry> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 3 {
        return Err(GdError::Parse {
            line: line_no,
            msg: "need type, id, parent_ids".into(),
        });
    }
    let kind: GdKind = cols[0]
        .parse()
        .map_err(|msg: String| GdError::Parse { line: line_no, msg })?;
    let id: u32 = cols[1].parse().map_err(|_| GdError::Parse {
        line: line_no,
        msg: "invalid id".into(),
    })?;
    let parent_ids = parse_parents(cols[2], line_no)?;
    let n_req = kind.field_count();
    let rest = &cols[3..];
    if rest.len() < n_req {
        return Err(GdError::Parse {
            line: line_no,
            msg: format!("{} needs {n_req} positional fields", kind.as_str()),
        });
    }
    let fields: Vec<String> = rest[..n_req].iter().map(|s| (*s).to_string()).collect();
    let mut attrs = BTreeMap::new();
    for extra in &rest[n_req..] {
        if let Some((k, v)) = extra.split_once('=') {
            attrs.insert(k.to_string(), v.to_string());
        } else if !extra.is_empty() {
            return Err(GdError::Parse {
                line: line_no,
                msg: format!("expected key=value, got {extra}"),
            });
        }
    }
    Ok(GdEntry {
        kind,
        id,
        parent_ids,
        fields,
        attrs,
    })
}

fn parse_parents(s: &str, line_no: usize) -> Result<Vec<u32>> {
    if s.is_empty() || s == "." {
        return Ok(Vec::new());
    }
    s.split(',')
        .map(|p| {
            p.parse().map_err(|_| GdError::Parse {
                line: line_no,
                msg: format!("invalid parent id {p}"),
            })
        })
        .collect()
}

impl fmt::Display for GenomeDiff {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_gd_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_gd() -> &'static str {
        "\
#=GENOME_DIFF\t1.0
#=AUTHOR\tprokdiff-test
SNP\t1\t.\tNC_000913\t100\tA
INS\t2\t.\tNC_000913\t200\tAT
DEL\t3\t.\tNC_000913\t300\t2
MOB\t4\t.\tNC_000913\t400\tIS1\t+\t8
JC\t5\t.\tNC_000913\t10\t+\tNC_000913\t500\t-\t0
UN\t6\t.\tNC_000913\t900\t910
"
    }

    #[test]
    fn parses_snp_ins_del_mob_jc() {
        let gd = GenomeDiff::parse(sample_gd()).expect("parse");
        assert_eq!(gd.entries.len(), 6);
        assert_eq!(gd.entries[0].kind, GdKind::Snp);
        assert_eq!(gd.entries[0].fields, ["NC_000913", "100", "A"]);
        assert_eq!(gd.entries[1].kind, GdKind::Ins);
        assert_eq!(gd.entries[1].fields[2], "AT");
        assert_eq!(gd.entries[2].kind, GdKind::Del);
        assert_eq!(gd.entries[2].fields[2], "2");
        assert_eq!(gd.entries[3].kind, GdKind::Mob);
        assert_eq!(gd.entries[3].fields[2], "IS1");
        assert_eq!(gd.entries[4].kind, GdKind::Jc);
        assert_eq!(gd.entries[4].fields[6], "0");
        assert_eq!(gd.entries[5].kind, GdKind::Un);
    }

    #[test]
    fn roundtrip_preserves_mutations() {
        let gd = GenomeDiff::parse(sample_gd()).unwrap();
        let again = GenomeDiff::parse(&gd.to_gd_string()).unwrap();
        assert_eq!(gd.entries, again.entries);
    }

    #[test]
    fn subtract_removes_matching_snp_allele() {
        let edited = GenomeDiff {
            metadata: vec![("GENOME_DIFF".into(), "1.0".into())],
            entries: vec![
                GdEntry::snp(1, "chr", 100, "A"),
                GdEntry::snp(2, "chr", 200, "T"),
            ],
        };
        let starter = GenomeDiff {
            metadata: vec![("GENOME_DIFF".into(), "1.0".into())],
            entries: vec![GdEntry::snp(9, "chr", 100, "A")],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].fields[1], "200");
    }

    #[test]
    fn subtract_keeps_snp_when_allele_differs() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::snp(1, "chr", 100, "A")],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::snp(1, "chr", 100, "C")],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
    }

    #[test]
    fn subtract_matches_ins_and_del_by_allele_or_size() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![
                GdEntry::ins(1, "chr", 50, "GG"),
                GdEntry::del(2, "chr", 80, 2),
                GdEntry::ins(3, "chr", 50, "TT"),
            ],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![
                GdEntry::ins(1, "chr", 50, "GG"),
                GdEntry::del(2, "chr", 80, 2),
            ],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].fields[2], "TT");
    }

    #[test]
    fn subtract_jc_requires_exact_coordinates() {
        // First-period product subtract matches gdtools SUBTRACT (exact fields).
        // A 1 bp junction jitter therefore does *not* cancel — documented risk.
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 100, "+", "chr", 500, "-", 0)],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 101, "+", "chr", 500, "-", 0)],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].fields[1], "100");
    }

    #[test]
    fn subtract_does_not_drop_un_as_oracle_failure() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![
                GdEntry {
                    kind: GdKind::Un,
                    id: 1,
                    parent_ids: vec![],
                    fields: vec!["chr".into(), "1".into(), "10".into()],
                    attrs: BTreeMap::new(),
                },
                GdEntry::snp(2, "chr", 20, "A"),
            ],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry {
                kind: GdKind::Un,
                id: 1,
                parent_ids: vec![],
                fields: vec!["chr".into(), "1".into(), "10".into()],
                attrs: BTreeMap::new(),
            }],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 2);
    }
}
