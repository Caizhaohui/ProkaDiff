//! Genome Diff subset I/O and `gdtools SUBTRACT`-style set difference.
//!
//! Matching key: mutation type + coordinates + allele (see `docs/schema.md`).

#![deny(unsafe_code)]

use std::cmp::Ordering;
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

pub mod normalize;
pub use normalize::{right_align_del, right_align_ins};

/// Three-letter mutation and two-letter evidence types in the first-period subset.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum GdKind {
    Snp,
    Sub,
    Ins,
    Del,
    Mob,
    Amp,
    Con,
    Inv,
    Jc,
    Un,
    Ra,
    Mc,
}

impl GdKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Snp => "SNP",
            Self::Sub => "SUB",
            Self::Ins => "INS",
            Self::Del => "DEL",
            Self::Mob => "MOB",
            Self::Amp => "AMP",
            Self::Con => "CON",
            Self::Inv => "INV",
            Self::Jc => "JC",
            Self::Un => "UN",
            Self::Ra => "RA",
            Self::Mc => "MC",
        }
    }

    fn field_count(self) -> usize {
        match self {
            Self::Snp | Self::Ins => 3,
            Self::Sub => 4,
            Self::Del | Self::Un | Self::Inv => 3,
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
            "SUB" => Self::Sub,
            "INS" => Self::Ins,
            "DEL" => Self::Del,
            "MOB" => Self::Mob,
            "AMP" => Self::Amp,
            "CON" => Self::Con,
            "INV" => Self::Inv,
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

    pub fn sub(
        id: u32,
        seq_id: impl Into<String>,
        position: u64,
        size: u64,
        new_seq: impl Into<String>,
    ) -> Self {
        Self {
            kind: GdKind::Sub,
            id,
            parent_ids: Vec::new(),
            fields: vec![
                seq_id.into(),
                position.to_string(),
                size.to_string(),
                new_seq.into(),
            ],
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

    pub fn ra(
        id: u32,
        seq_id: impl Into<String>,
        position: u64,
        insert_position: u64,
        ref_base: impl Into<String>,
        new_base: impl Into<String>,
    ) -> Self {
        Self {
            kind: GdKind::Ra,
            id,
            parent_ids: Vec::new(),
            fields: vec![
                seq_id.into(),
                position.to_string(),
                insert_position.to_string(),
                ref_base.into(),
                new_base.into(),
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

    pub fn amp(
        id: u32,
        seq_id: impl Into<String>,
        position: u64,
        size: u64,
        new_copy_number: u64,
    ) -> Self {
        Self {
            kind: GdKind::Amp,
            id,
            parent_ids: Vec::new(),
            fields: vec![
                seq_id.into(),
                position.to_string(),
                size.to_string(),
                new_copy_number.to_string(),
            ],
            attrs: BTreeMap::new(),
        }
    }

    pub fn con(
        id: u32,
        seq_id: impl Into<String>,
        position: u64,
        size: u64,
        region: impl Into<String>,
    ) -> Self {
        Self {
            kind: GdKind::Con,
            id,
            parent_ids: Vec::new(),
            fields: vec![
                seq_id.into(),
                position.to_string(),
                size.to_string(),
                region.into(),
            ],
            attrs: BTreeMap::new(),
        }
    }

    pub fn inv(id: u32, seq_id: impl Into<String>, position: u64, size: u64) -> Self {
        Self {
            kind: GdKind::Inv,
            id,
            parent_ids: Vec::new(),
            fields: vec![seq_id.into(), position.to_string(), size.to_string()],
            attrs: BTreeMap::new(),
        }
    }

    pub fn new_inv(id: u32, seq_id: impl Into<String>, position: u64, size: u64) -> Self {
        Self::inv(id, seq_id, position, size)
    }

    pub fn del_size(&self) -> Option<u64> {
        if self.kind == GdKind::Del {
            self.fields.get(2).and_then(|s| s.parse().ok())
        } else {
            None
        }
    }

    pub fn inv_size(&self) -> Option<u64> {
        if self.kind == GdKind::Inv {
            self.fields.get(2).and_then(|s| s.parse().ok())
        } else {
            None
        }
    }

    pub fn sub_size(&self) -> Option<u64> {
        if self.kind == GdKind::Sub {
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

pub const DEFAULT_JC_SUBTRACT_TOL_BP: u64 = 5;
pub const DEFAULT_MOB_SUBTRACT_TOL_BP: u64 = 5;
pub const DEFAULT_DEL_SUBTRACT_TOL_BP: u64 = 5;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CanonicalJcSide {
    pub seq_id: String,
    pub strand: String,
    pub pos: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CanonicalJc {
    pub side1: CanonicalJcSide,
    pub side2: CanonicalJcSide,
    pub overlap: i64,
}

fn normalize_strand(s: &str) -> String {
    match s {
        "+" | "1" => "+".to_string(),
        "-" | "-1" => "-".to_string(),
        other => other.to_string(),
    }
}

impl CanonicalJc {
    pub fn parse(fields: &[String]) -> Option<Self> {
        if fields.len() < 6 {
            return None;
        }
        let pos1: u64 = fields[1].parse().ok()?;
        let pos2: u64 = fields[4].parse().ok()?;
        let s1 = CanonicalJcSide {
            seq_id: fields[0].clone(),
            strand: normalize_strand(&fields[2]),
            pos: pos1,
        };
        let s2 = CanonicalJcSide {
            seq_id: fields[3].clone(),
            strand: normalize_strand(&fields[5]),
            pos: pos2,
        };
        let (side1, side2) = if s1 <= s2 { (s1, s2) } else { (s2, s1) };
        let overlap = fields.get(6).and_then(|s| s.parse().ok())?;
        Some(Self {
            side1,
            side2,
            overlap,
        })
    }

    pub fn matches_tolerant(&self, other: &Self, tol: u64) -> bool {
        self.overlap == other.overlap
            && self.side1.seq_id == other.side1.seq_id
            && self.side1.strand == other.side1.strand
            && self.side1.pos.abs_diff(other.side1.pos) <= tol
            && self.side2.seq_id == other.side2.seq_id
            && self.side2.strand == other.side2.strand
            && self.side2.pos.abs_diff(other.side2.pos) <= tol
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CanonicalMob {
    pub seq_id: String,
    pub pos: u64,
    pub repeat_name: String,
    pub strand: String,
    pub duplication_size: i64,
}

impl CanonicalMob {
    pub fn parse(fields: &[String]) -> Option<Self> {
        if fields.len() < 5 {
            return None;
        }
        let pos: u64 = fields[1].parse().ok()?;
        let duplication_size: i64 = fields[4].parse().ok()?;
        Some(Self {
            seq_id: fields[0].clone(),
            pos,
            repeat_name: fields[2].clone(),
            strand: normalize_strand(&fields[3]),
            duplication_size,
        })
    }

    pub fn matches_tolerant(&self, other: &Self, tol: u64) -> bool {
        self.seq_id == other.seq_id
            && self.repeat_name == other.repeat_name
            && self.strand == other.strand
            && self.duplication_size == other.duplication_size
            && self.pos.abs_diff(other.pos) <= tol
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CanonicalDel {
    pub seq_id: String,
    pub start: u64,
    pub size: u64,
}

impl CanonicalDel {
    pub fn parse(fields: &[String]) -> Option<Self> {
        if fields.len() < 3 {
            return None;
        }
        let start: u64 = fields[1].parse().ok()?;
        let size: u64 = fields[2].parse().ok()?;
        Some(Self {
            seq_id: fields[0].clone(),
            start,
            size,
        })
    }

    pub fn end(&self) -> Option<u64> {
        self.start.checked_add(self.size.checked_sub(1)?)
    }

    pub fn matches_tolerant(&self, other: &Self, tol: u64) -> bool {
        if self.seq_id != other.seq_id {
            return false;
        }
        if self.size <= 2 || other.size <= 2 {
            return self.start == other.start && self.size == other.size;
        }
        let Some(self_end) = self.end() else {
            return false;
        };
        let Some(other_end) = other.end() else {
            return false;
        };
        self.start.abs_diff(other.start) <= tol && self_end.abs_diff(other_end) <= tol
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct MatchCost {
    max_delta: u64,
    sum_delta: u64,
}

impl MatchCost {
    const ZERO: Self = Self {
        max_delta: 0,
        sum_delta: 0,
    };

    fn checked_add(self, other: Self) -> Option<Self> {
        Some(Self {
            max_delta: self.max_delta.checked_add(other.max_delta)?,
            sum_delta: self.sum_delta.checked_add(other.sum_delta)?,
        })
    }
}

impl Ord for MatchCost {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.max_delta, self.sum_delta).cmp(&(other.max_delta, other.sum_delta))
    }
}

impl PartialOrd for MatchCost {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MatchCandidate {
    edited_index: usize,
    starter_index: usize,
    cost: MatchCost,
    edited_key: Vec<u8>,
    starter_key: Vec<u8>,
}

impl Ord for MatchCandidate {
    fn cmp(&self, other: &Self) -> Ordering {
        (
            self.cost,
            self.edited_key.as_slice(),
            self.starter_key.as_slice(),
            self.edited_index,
            self.starter_index,
        )
            .cmp(&(
                other.cost,
                other.edited_key.as_slice(),
                other.starter_key.as_slice(),
                other.edited_index,
                other.starter_index,
            ))
    }
}

impl PartialOrd for MatchCandidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Copy, Debug)]
struct ResidualEdge {
    to: usize,
    reverse_index: usize,
    capacity: u8,
    cost: MatchCost,
}

fn add_match_edge(graph: &mut [Vec<ResidualEdge>], from: usize, to: usize, cost: MatchCost) {
    let forward_index = graph[from].len();
    let reverse_index = graph[to].len();
    graph[from].push(ResidualEdge {
        to,
        reverse_index,
        capacity: 1,
        cost,
    });
    graph[to].push(ResidualEdge {
        to: from,
        reverse_index: forward_index,
        capacity: 0,
        cost,
    });
}

fn minimum_cost_maximum_matching(
    edited_count: usize,
    starter_count: usize,
    candidates: &[MatchCandidate],
) -> Vec<(usize, usize)> {
    let source = 0;
    let edited_start = 1;
    let starter_start = edited_start + edited_count;
    let sink = starter_start + starter_count;
    let mut graph = vec![Vec::new(); sink + 1];

    for edited_index in 0..edited_count {
        add_match_edge(
            &mut graph,
            source,
            edited_start + edited_index,
            MatchCost::ZERO,
        );
    }
    for candidate in candidates {
        add_match_edge(
            &mut graph,
            edited_start + candidate.edited_index,
            starter_start + candidate.starter_index,
            candidate.cost,
        );
    }
    for starter_index in 0..starter_count {
        add_match_edge(
            &mut graph,
            starter_start + starter_index,
            sink,
            MatchCost::ZERO,
        );
    }

    loop {
        let mut distance = vec![None; graph.len()];
        let mut previous = vec![None; graph.len()];
        distance[source] = Some(MatchCost::ZERO);

        for _ in 1..graph.len() {
            let mut changed = false;
            for node in 0..graph.len() {
                let Some(prefix_cost) = distance[node] else {
                    continue;
                };
                for (edge_index, edge) in graph[node].iter().enumerate() {
                    if edge.capacity == 0 {
                        continue;
                    }
                    let Some(candidate_cost) = prefix_cost.checked_add(edge.cost) else {
                        continue;
                    };
                    if distance[edge.to].is_none_or(|current| candidate_cost < current) {
                        distance[edge.to] = Some(candidate_cost);
                        previous[edge.to] = Some((node, edge_index));
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }

        if distance[sink].is_none() {
            break;
        }

        let mut node = sink;
        while node != source {
            let Some((from, edge_index)) = previous[node] else {
                return collect_min_cost_matches(
                    &graph,
                    edited_start,
                    starter_start,
                    starter_count,
                );
            };
            let reverse_index = graph[from][edge_index].reverse_index;
            graph[from][edge_index].capacity = 0;
            graph[node][reverse_index].capacity = 1;
            node = from;
        }
    }

    collect_min_cost_matches(&graph, edited_start, starter_start, starter_count)
}

fn collect_min_cost_matches(
    graph: &[Vec<ResidualEdge>],
    edited_start: usize,
    starter_start: usize,
    starter_count: usize,
) -> Vec<(usize, usize)> {
    let starter_end = starter_start + starter_count;
    let mut matches = Vec::new();
    for (edited_index, edges) in graph[edited_start..starter_start].iter().enumerate() {
        for edge in edges {
            if (starter_start..starter_end).contains(&edge.to) && edge.capacity == 0 {
                matches.push((edited_index, edge.to - starter_start));
            }
        }
    }
    matches
}

fn canonical_subtract_key(entry: &GdEntry) -> Option<Vec<u8>> {
    match entry.kind {
        GdKind::Un => None,
        GdKind::Jc => {
            let jc = CanonicalJc::parse(&entry.fields)?;
            Some(
                format!(
                    "JC|{}|{}|{}|{}|{}|{}|{}",
                    jc.side1.seq_id,
                    jc.side1.strand,
                    jc.side1.pos,
                    jc.side2.seq_id,
                    jc.side2.strand,
                    jc.side2.pos,
                    jc.overlap
                )
                .into_bytes(),
            )
        }
        GdKind::Mob => {
            let mob = CanonicalMob::parse(&entry.fields)?;
            Some(
                format!(
                    "MOB|{}|{}|{}|{}|{}",
                    mob.seq_id, mob.pos, mob.repeat_name, mob.strand, mob.duplication_size
                )
                .into_bytes(),
            )
        }
        _ => Some(entry.subtract_key().into_bytes()),
    }
}

fn tolerant_candidate_cost(
    edited: &GdEntry,
    starter: &GdEntry,
    tol: Tolerances,
) -> Option<MatchCost> {
    let deltas = match (edited.kind, starter.kind) {
        (GdKind::Jc, GdKind::Jc) => {
            let edited_jc = CanonicalJc::parse(&edited.fields)?;
            let starter_jc = CanonicalJc::parse(&starter.fields)?;
            if !edited_jc.matches_tolerant(&starter_jc, tol.jc) {
                return None;
            }
            [
                edited_jc.side1.pos.abs_diff(starter_jc.side1.pos),
                edited_jc.side2.pos.abs_diff(starter_jc.side2.pos),
            ]
        }
        (GdKind::Mob, GdKind::Mob) => {
            let edited_mob = CanonicalMob::parse(&edited.fields)?;
            let starter_mob = CanonicalMob::parse(&starter.fields)?;
            if !edited_mob.matches_tolerant(&starter_mob, tol.mob) {
                return None;
            }
            let delta = edited_mob.pos.abs_diff(starter_mob.pos);
            [delta, delta]
        }
        (GdKind::Del, GdKind::Del) => {
            let edited_del = CanonicalDel::parse(&edited.fields)?;
            let starter_del = CanonicalDel::parse(&starter.fields)?;
            if edited_del.size <= 2 || starter_del.size <= 2 {
                return None;
            }
            if !edited_del.matches_tolerant(&starter_del, tol.del) {
                return None;
            }
            [
                edited_del.start.abs_diff(starter_del.start),
                edited_del.end()?.abs_diff(starter_del.end()?),
            ]
        }
        (GdKind::Snp, _)
        | (GdKind::Sub, _)
        | (GdKind::Ins, _)
        | (GdKind::Amp, _)
        | (GdKind::Con, _)
        | (GdKind::Inv, _)
        | (GdKind::Un, _)
        | (GdKind::Ra, _)
        | (GdKind::Mc, _) => return None,
        (GdKind::Del, _) | (GdKind::Mob, _) | (GdKind::Jc, _) => return None,
    };
    Some(MatchCost {
        max_delta: deltas[0].max(deltas[1]),
        sum_delta: deltas[0].checked_add(deltas[1])?,
    })
}

#[derive(Clone, Copy, Debug)]
struct Tolerances {
    jc: u64,
    mob: u64,
    del: u64,
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
}

impl FromStr for GenomeDiff {
    type Err = GdError;

    fn from_str(s: &str) -> Result<Self> {
        Self::parse(s)
    }
}

impl GenomeDiff {
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
    /// oracle failure). Exact string match on kind and fields.
    pub fn subtract_exact(&self, other: &GenomeDiff) -> GenomeDiff {
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

    /// Tolerant set subtraction: applies coordinate tolerance to JC, MOB, and structural DEL
    /// entries to avoid false-positive unintended calls when independent runs
    /// pick representative coordinates that jitter within clustering tolerance (±5 bp).
    /// Non-JC/MOB/DEL entries use exact subtraction keys.
    pub fn subtract_tolerant(
        &self,
        other: &GenomeDiff,
        jc_tol_bp: u64,
        mob_tol_bp: u64,
        del_tol_bp: u64,
    ) -> GenomeDiff {
        let tolerances = Tolerances {
            jc: jc_tol_bp.min(DEFAULT_JC_SUBTRACT_TOL_BP),
            mob: mob_tol_bp.min(DEFAULT_MOB_SUBTRACT_TOL_BP),
            del: del_tol_bp.min(DEFAULT_DEL_SUBTRACT_TOL_BP),
        };

        let mut used_other = vec![false; other.entries.len()];
        let mut removed_self = vec![false; self.entries.len()];

        let mut edited_exact: Vec<(Vec<u8>, usize)> = self
            .entries
            .iter()
            .enumerate()
            .filter_map(|(idx, entry)| canonical_subtract_key(entry).map(|key| (key, idx)))
            .collect();
        let mut starter_exact: Vec<(Vec<u8>, usize)> = other
            .entries
            .iter()
            .enumerate()
            .filter_map(|(idx, entry)| canonical_subtract_key(entry).map(|key| (key, idx)))
            .collect();
        edited_exact.sort_unstable();
        starter_exact.sort_unstable();

        let mut starter_cursor = 0usize;
        for (edited_key, edited_index) in edited_exact {
            while starter_cursor < starter_exact.len()
                && starter_exact[starter_cursor].0 < edited_key
            {
                starter_cursor += 1;
            }
            let mut probe = starter_cursor;
            while probe < starter_exact.len() && starter_exact[probe].0 == edited_key {
                let starter_index = starter_exact[probe].1;
                if !used_other[starter_index] {
                    used_other[starter_index] = true;
                    removed_self[edited_index] = true;
                    break;
                }
                probe += 1;
            }
        }

        let mut candidates = Vec::new();
        for (edited_index, edited) in self.entries.iter().enumerate() {
            if removed_self[edited_index] {
                continue;
            }
            let Some(edited_key) = canonical_subtract_key(edited) else {
                continue;
            };
            for (starter_index, starter) in other.entries.iter().enumerate() {
                if used_other[starter_index] {
                    continue;
                }
                let Some(starter_key) = canonical_subtract_key(starter) else {
                    continue;
                };
                if let Some(cost) = tolerant_candidate_cost(edited, starter, tolerances) {
                    candidates.push(MatchCandidate {
                        edited_index,
                        starter_index,
                        cost,
                        edited_key: edited_key.clone(),
                        starter_key,
                    });
                }
            }
        }
        candidates.sort_unstable();
        for (edited_index, starter_index) in
            minimum_cost_maximum_matching(self.entries.len(), other.entries.len(), &candidates)
        {
            if !removed_self[edited_index] && !used_other[starter_index] {
                removed_self[edited_index] = true;
                used_other[starter_index] = true;
            }
        }

        let entries = self
            .entries
            .iter()
            .enumerate()
            .filter(|(idx, entry)| entry.kind == GdKind::Un || !removed_self[*idx])
            .map(|(_, entry)| entry.clone())
            .collect();

        GenomeDiff {
            metadata: self.metadata.clone(),
            entries,
        }
    }

    /// ProkaDiff product subtract: cancels matching mutations, with ±5 bp
    /// tolerance for JC, MOB, and structural DEL coordinates to accommodate independent clustering jitter.
    pub fn subtract(&self, other: &GenomeDiff) -> GenomeDiff {
        self.subtract_tolerant(
            other,
            DEFAULT_JC_SUBTRACT_TOL_BP,
            DEFAULT_MOB_SUBTRACT_TOL_BP,
            DEFAULT_DEL_SUBTRACT_TOL_BP,
        )
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
    fn subtract_exact_requires_exact_coordinates() {
        // Exact subtract matches gdtools SUBTRACT (exact fields).
        // A 1 bp junction jitter therefore does *not* cancel under subtract_exact.
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 100, "+", "chr", 500, "-", 0)],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 101, "+", "chr", 500, "-", 0)],
        };
        let out = edited.subtract_exact(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].fields[1], "100");
    }

    #[test]
    fn subtract_tolerant_cancels_jittered_jc_within_tolerance() {
        // Product subtract cancels jittered JCs within default 5 bp tolerance
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![
                GdEntry::jc(1, "chr", 100, "+", "chr", 500, "-", 0),
                GdEntry::jc(2, "chr", 2000, "+", "chr", 3000, "-", 0),
            ],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 102, "+", "chr", 501, "-", 0)],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].fields[1], "2000");
    }

    #[test]
    fn subtract_tolerant_preserves_distinct_jc_beyond_tolerance() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 100, "+", "chr", 500, "-", 0)],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 108, "+", "chr", 500, "-", 0)],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].fields[1], "100");
    }

    #[test]
    fn subtract_tolerant_cancels_jittered_mob() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![
                GdEntry::mob(1, "chr", 601, "IS150", "+", 3),
                GdEntry::snp(2, "chr", 800, "A"),
            ],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::mob(1, "chr", 603, "IS150", "+", 3)],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].kind, GdKind::Snp);
    }

    #[test]
    fn subtract_tolerant_cancels_jittered_structural_del() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![
                GdEntry::del(1, "chr", 4999, 304),
                GdEntry::snp(2, "chr", 8000, "A"),
            ],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::del(1, "chr", 5000, 300)],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].kind, GdKind::Snp);
    }

    #[test]
    fn subtract_tolerant_preserves_distinct_del_beyond_tolerance() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::del(1, "chr", 5000, 300)],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::del(1, "chr", 5020, 300)],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
    }

    #[test]
    fn subtract_tolerant_handles_inverted_side_order() {
        // side1 and side2 order may be flipped between independent callers
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 500, "-", "chr", 100, "+", 0)],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::jc(1, "chr", 101, "+", "chr", 502, "-", 0)],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 0);
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

    #[test]
    fn parses_sub_and_roundtrip() {
        let line = "SUB\t1\t2,3\tchr\t201\t2\tTT\tfoo=bar\n";
        let gd = GenomeDiff::parse(line).expect("should parse SUB");
        assert_eq!(gd.entries.len(), 1);
        let e = &gd.entries[0];
        assert_eq!(e.kind, GdKind::Sub);
        assert_eq!(e.id, 1);
        assert_eq!(e.parent_ids, vec![2, 3]);
        assert_eq!(e.seq_id(), Some("chr"));
        assert_eq!(e.position(), Some(201));
        assert_eq!(e.sub_size(), Some(2));
        assert_eq!(e.fields.get(3).map(String::as_str), Some("TT"));
        assert_eq!(e.attrs.get("foo").map(String::as_str), Some("bar"));
        assert_eq!(e.to_line(), "SUB\t1\t2,3\tchr\t201\t2\tTT\tfoo=bar");
    }

    #[test]
    fn subtract_removes_matching_sub() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![
                GdEntry::sub(1, "chr", 100, 2, "TT"),
                GdEntry::snp(2, "chr", 300, "G"),
            ],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::sub(1, "chr", 100, 2, "TT")],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].kind, GdKind::Snp);
    }

    #[test]
    fn subtract_keeps_sub_when_allele_differs() {
        let edited = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::sub(1, "chr", 100, 2, "TT")],
        };
        let starter = GenomeDiff {
            metadata: vec![],
            entries: vec![GdEntry::sub(1, "chr", 100, 2, "TG")],
        };
        let out = edited.subtract(&starter);
        assert_eq!(out.entries.len(), 1);
        assert_eq!(out.entries[0].kind, GdKind::Sub);
    }

    #[test]
    fn test_inv_roundtrip() {
        let text = "#=GENOME_DIFF\t1.0\nINV\t10\t.\tchr1\t1000\t500\n";
        let gd = GenomeDiff::from_str(text).unwrap();
        assert_eq!(gd.entries.len(), 1);
        let e = &gd.entries[0];
        assert_eq!(e.kind, GdKind::Inv);
        assert_eq!(e.seq_id(), Some("chr1"));
        assert_eq!(e.position(), Some(1000));
        assert_eq!(e.inv_size(), Some(500));
        assert_eq!(gd.to_string(), text);
    }

    #[test]
    fn inv_rejects_missing_positional_fields() {
        let text = "INV\t10\t.\tchr1\t1000\n";
        assert!(GenomeDiff::from_str(text).is_err());
    }

    #[test]
    fn inv_rejects_invalid_type() {
        let text = "INVX\t10\t.\tchr1\t1000\t500\n";
        assert!(GenomeDiff::from_str(text).is_err());
    }
}
