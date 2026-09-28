//! Internal repeat-family junction evidence.
//!
//! These records preserve observed junction support when a repeated reference
//! sequence prevents assignment to one copy. They are diagnostic evidence only:
//! unresolved records are never converted into a GenomeDiff mutation.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt::Write as _;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::fasta::FastaRecord;
use crate::jc::{accept_junction, JunctionSupport};

/// Safe diagnostic retention limit. Above this, the exact clip sequence plus
/// placement count remains a lossless reference query, but individual copies
/// are not materialized and the result cannot resolve to an exact junction.
pub const MAX_RETAINED_REPEAT_PLACEMENTS: usize = 4096;
/// A reciprocal search may inspect this many indexed families for one seed.
/// Crossing the bound is diagnostic-only: resolution is never attempted from
/// a prefix of the candidates.
const MAX_RECIPROCAL_CANDIDATES_PER_SEED: usize = 5_000;
/// Bound the reverse index itself so pathological repeat families remain
/// observable without permitting unbounded diagnostic work.
const MAX_RECIPROCAL_INDEX_ENTRIES: usize = 200_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum RepeatAmbiguousState {
    Candidate,
    SupportedAmbiguous,
    Resolved,
    ResourceLimit,
}

impl RepeatAmbiguousState {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Candidate => "CANDIDATE",
            Self::SupportedAmbiguous => "SUPPORTED_AMBIGUOUS",
            Self::Resolved => "RESOLVED",
            Self::ResourceLimit => "RESOURCE_LIMIT",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum RepeatEvidenceSource {
    PrimaryClip,
    Stage2Mosaic,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum CopyResolutionSource {
    UniqueFlankExtension,
    UniqueMateGeometry,
    ReciprocalBreakpoint,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CopyResolutionEvidence {
    pub source: CopyResolutionSource,
    pub molecule_id: u64,
    pub copy: RepeatCopyPlacement,
}

pub const MATE_UNKNOWN: u8 = 0;
pub const MATE_1: u8 = 1;
pub const MATE_2: u8 = 2;

/// Canonical mate derivation from SAM flags and read name.
///
/// Precedence is explicit and consistent across all parsing paths:
/// 1. SAM flags when segmented:
///    - first segment -> 1 (MATE_1)
///    - last segment -> 2 (MATE_2)
/// 2. Explicit terminal suffix in read name:
///    - terminal "/1" -> 1 (MATE_1)
///    - terminal "/2" -> 2 (MATE_2)
/// 3. Otherwise -> 0 (MATE_UNKNOWN)
///
/// Caller-specific context/file roles never invent mate 1/2.
pub fn canonical_mate(
    flags: Option<&noodles::sam::alignment::record::Flags>,
    raw_name: &str,
) -> u8 {
    if let Some(flags) = flags {
        if flags.is_first_segment() {
            return MATE_1;
        } else if flags.is_last_segment() {
            return MATE_2;
        }
    }
    if raw_name.ends_with("/1") {
        MATE_1
    } else if raw_name.ends_with("/2") {
        MATE_2
    } else {
        MATE_UNKNOWN
    }
}

/// Canonical read-molecule identity shared by primary clipping, stage-2
/// alignment, and reciprocal evidence. `qname` has already had a terminal
/// `/1` or `/2` normalised by the alignment reader; mate remains explicit so
/// paired reads cannot collapse into one molecule.
pub fn normalized_molecule_id(qname: &str, mate: u8) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    qname.hash(&mut hasher);
    mate.hash(&mut hasher);
    hasher.finish()
}

/// Canonical molecule identity constructor.
///
/// Normalizes:
/// - QNAME (stripping `/1` or `/2`)
/// - mate code (via [`canonical_mate`])
///
/// Returns (normalized_qname, canonical_mate, canonical_molecule_id).
pub fn canonical_molecule_identity(
    raw_name: &str,
    flags: Option<&noodles::sam::alignment::record::Flags>,
) -> (String, u8, u64) {
    let mate = canonical_mate(flags, raw_name);
    let norm_qname = crate::align::normalize_qname(raw_name).to_string();
    let id = normalized_molecule_id(&norm_qname, mate);
    (norm_qname, mate, id)
}

impl RepeatEvidenceSource {
    const fn as_str(self) -> &'static str {
        match self {
            Self::PrimaryClip => "PRIMARY_CLIP",
            Self::Stage2Mosaic => "STAGE2_MOSAIC",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct RepeatCopyPlacement {
    pub contig_idx: usize,
    pub position_1: u64,
    pub minus: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepeatAmbiguousSeed {
    pub source: RepeatEvidenceSource,
    pub anchor_contig_idx: usize,
    pub anchor_position_1: u64,
    pub anchor_minus: bool,
    /// Raw observed placements (including trivial continuations) for diagnostic provenance.
    pub observed_placements: Vec<RepeatCopyPlacement>,
    /// Scientifically eligible feasible copies (excluding trivial continuations and overflowed subsets).
    pub eligible_copies: Vec<RepeatCopyPlacement>,
    pub placement_count: usize,
    pub placement_family: Vec<u8>,
    pub placement_sequence: Vec<u8>,
    pub placement_set_complete: bool,
    pub resource_complete: bool,
    pub resolution_complete: bool,
    pub overlap: i64,
    pub molecule_id: u64,
    pub molecule_minus: bool,
    pub clip_length: usize,
    pub aligned_length: usize,
    pub effective_mapq: u8,
    pub alignment_score: Option<i32>,
    pub pair_geometry: &'static str,
    pub unique_anchor_qualified: bool,
    pub copy_resolved: bool,
    pub resolution_molecule_id: Option<u64>,
    pub reciprocal_evidence: bool,
    pub rejection_reason: Option<&'static str>,
    /// A resource bound stopped family analysis. This is distinct from an
    /// incomplete placement set: the observed evidence remains complete as
    /// far as placement search is concerned, but it cannot resolve exactly.
    pub resource_limit_reason: Option<&'static str>,
}

pub fn resolve_repeat_copy(
    seed: &RepeatAmbiguousSeed,
    evidence: &[CopyResolutionEvidence],
) -> RepeatAmbiguousSeed {
    if !seed.placement_set_complete
        || !seed.resource_complete
        || seed.resource_limit_reason.is_some()
        || seed.eligible_copies.len() < 2
    {
        return seed.clone();
    }
    let feasible: BTreeSet<RepeatCopyPlacement> = seed.eligible_copies.iter().copied().collect();
    let eligible: Vec<CopyResolutionEvidence> = evidence
        .iter()
        .copied()
        .filter(|item| item.molecule_id != seed.molecule_id && feasible.contains(&item.copy))
        .collect();
    let targets: BTreeSet<RepeatCopyPlacement> = eligible.iter().map(|item| item.copy).collect();
    if targets.len() != 1 {
        return seed.clone();
    }
    let Some(copy) = targets.into_iter().next() else {
        return seed.clone();
    };
    let Some(resolution) = eligible.first() else {
        return seed.clone();
    };
    let mut resolved = seed.clone();
    resolved.eligible_copies = vec![copy];
    resolved.copy_resolved = true;
    resolved.resolution_complete = true;
    resolved.resolution_molecule_id = Some(resolution.molecule_id);
    resolved.reciprocal_evidence = eligible
        .iter()
        .any(|item| item.source == CopyResolutionSource::ReciprocalBreakpoint);
    resolved
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepeatAmbiguousJunctionEvidence {
    pub state: RepeatAmbiguousState,
    pub reference_contigs: Vec<String>,
    pub anchor_contig_idx: usize,
    pub anchor_position_1: u64,
    pub anchor_minus: bool,
    pub feasible_copies: Vec<RepeatCopyPlacement>,
    pub observed_placements: Vec<RepeatCopyPlacement>,
    pub placement_count: usize,
    pub placement_family: Vec<u8>,
    pub placement_sequence: Vec<u8>,
    pub placement_set_complete: bool,
    pub resource_complete: bool,
    pub resolution_complete: bool,
    pub breakpoint_start_1: u64,
    pub breakpoint_end_1: u64,
    pub overlap: i64,
    pub molecule_ids: Vec<u64>,
    pub plus_molecules: usize,
    pub minus_molecules: usize,
    pub clip_lengths: Vec<usize>,
    pub aligned_lengths: Vec<usize>,
    pub effective_mapqs: Vec<u8>,
    pub alignment_scores: Vec<i32>,
    pub sources: Vec<RepeatEvidenceSource>,
    pub pair_geometries: Vec<&'static str>,
    pub reciprocal_evidence: bool,
    pub resolution_molecule_ids: Vec<u64>,
    pub coverage_context: &'static str,
    pub rejection_reasons: Vec<&'static str>,
    pub exact_breakpoint_supported: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RepeatAmbiguousDiagnostics {
    pub records: Vec<RepeatAmbiguousJunctionEvidence>,
    pub candidate_count: usize,
    pub repeat_family_count: usize,
    pub overflow_count: usize,
    pub discard_count: usize,
    pub trivial_continuation_count: usize,
}

type FamilyKey = (
    usize,
    u64,
    bool,
    i64,
    bool,
    Vec<RepeatCopyPlacement>,
    Vec<u8>,
);

fn family_key(seed: &RepeatAmbiguousSeed) -> FamilyKey {
    if seed.unique_anchor_qualified
        && seed.placement_set_complete
        && seed.resource_complete
        && seed.resource_limit_reason.is_none()
    {
        let copies = seed
            .eligible_copies
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        (
            seed.anchor_contig_idx,
            seed.anchor_position_1,
            seed.anchor_minus,
            seed.overlap,
            true,
            copies,
            Vec::new(),
        )
    } else {
        (
            seed.anchor_contig_idx,
            seed.anchor_position_1,
            seed.anchor_minus,
            seed.overlap,
            false,
            Vec::new(),
            seed.placement_family.clone(),
        )
    }
}

impl RepeatAmbiguousDiagnostics {
    pub fn from_seeds(fasta: &[FastaRecord], seeds: Vec<RepeatAmbiguousSeed>) -> Self {
        // Trivial continuations are scientifically rejected: do not promote to CANDIDATE.
        let mut trivial_continuation_count = 0usize;
        let mut active_seeds = Vec::with_capacity(seeds.len());
        for seed in seeds {
            if seed.rejection_reason == Some("TRIVIAL_CONTINUATION_ONLY") {
                trivial_continuation_count += 1;
            } else {
                active_seeds.push(seed);
            }
        }
        let mut seeds = active_seeds;

        // Resource overflow is FAMILY-SCOPED.
        // If any seed in a family has a resource limit or incomplete placement set,
        // the whole affected family must become RESOURCE_LIMIT / non-resolvable.
        let mut incomplete_families: HashSet<Vec<u8>> = HashSet::new();
        for seed in &seeds {
            if seed.resource_limit_reason.is_some()
                || !seed.placement_set_complete
                || !seed.resource_complete
            {
                incomplete_families.insert(seed.placement_family.clone());
            }
        }
        if !incomplete_families.is_empty() {
            for seed in &mut seeds {
                if incomplete_families.contains(&seed.placement_family) {
                    if seed.resource_limit_reason.is_none() {
                        seed.resource_limit_reason = Some("RESOURCE_LIMIT_INCOMPLETE_FAMILY");
                    }
                    seed.placement_set_complete = false;
                    seed.resource_complete = false;
                    seed.resolution_complete = false;
                    seed.eligible_copies.clear();
                }
            }
        }

        let seeds = resolve_reciprocal_seeds(seeds);
        let mut families: BTreeMap<FamilyKey, Vec<RepeatAmbiguousSeed>> = BTreeMap::new();
        for seed in seeds {
            families.entry(family_key(&seed)).or_default().push(seed);
        }

        let mut diagnostics = Self {
            candidate_count: families.len(),
            repeat_family_count: families.len(),
            trivial_continuation_count,
            ..Self::default()
        };
        diagnostics.records = families
            .into_values()
            .map(|family| Self::aggregate_family(fasta, family))
            .collect();
        diagnostics.records.sort_by_key(|record| {
            (
                record.anchor_contig_idx,
                record.anchor_position_1,
                record.anchor_minus,
                record.placement_family.clone(),
            )
        });
        diagnostics.overflow_count = diagnostics
            .records
            .iter()
            .filter(|record| record.state == RepeatAmbiguousState::ResourceLimit)
            .count();
        diagnostics.discard_count = diagnostics
            .records
            .iter()
            .filter(|record| record.state == RepeatAmbiguousState::Candidate)
            .count();
        diagnostics
    }

    fn aggregate_family(
        fasta: &[FastaRecord],
        mut family: Vec<RepeatAmbiguousSeed>,
    ) -> RepeatAmbiguousJunctionEvidence {
        family.sort_by_key(|seed| {
            (
                seed.molecule_id,
                seed.source,
                seed.molecule_minus,
                seed.clip_length,
                seed.aligned_length,
            )
        });
        let first = &family[0];
        let mut copies = BTreeSet::new();
        let mut observed = BTreeSet::new();
        let mut molecule_support: BTreeMap<u64, (bool, usize, usize)> = BTreeMap::new();
        let mut clip_lengths = BTreeSet::new();
        let mut aligned_lengths = BTreeSet::new();
        let mut effective_mapqs = BTreeSet::new();
        let mut scores = BTreeSet::new();
        let mut sources = BTreeSet::new();
        let mut geometry = BTreeSet::new();
        let mut rejection_reasons = BTreeSet::new();
        let mut placement_count = 0usize;
        let mut complete = true;
        let mut resource_complete = true;
        let mut resolution_complete = false;
        let mut anchors_qualified = true;
        let mut copy_resolved = false;
        let mut resolution_molecules = BTreeSet::new();
        let mut reciprocal_evidence = false;

        for seed in &family {
            placement_count = placement_count.max(seed.placement_count);
            complete &= seed.placement_set_complete;
            resource_complete &= seed.resource_complete;
            resolution_complete |= seed.resolution_complete;
            anchors_qualified &= seed.unique_anchor_qualified;
            copy_resolved |= seed.copy_resolved;
            if let Some(molecule_id) = seed.resolution_molecule_id {
                resolution_molecules.insert(molecule_id);
            }
            reciprocal_evidence |= seed.reciprocal_evidence;
            copies.extend(seed.eligible_copies.iter().copied());
            observed.extend(seed.observed_placements.iter().copied());
            molecule_support
                .entry(seed.molecule_id)
                .and_modify(|support| {
                    support.1 = support.1.max(seed.aligned_length);
                    support.2 = support.2.max(seed.clip_length);
                })
                .or_insert((seed.molecule_minus, seed.aligned_length, seed.clip_length));
            clip_lengths.insert(seed.clip_length);
            aligned_lengths.insert(seed.aligned_length);
            effective_mapqs.insert(seed.effective_mapq);
            if let Some(score) = seed.alignment_score {
                scores.insert(score);
            }
            sources.insert(seed.source);
            geometry.insert(seed.pair_geometry);
            if let Some(reason) = seed.rejection_reason {
                rejection_reasons.insert(reason);
            }
            if let Some(reason) = seed.resource_limit_reason {
                rejection_reasons.insert(reason);
            }
        }

        let mut support = JunctionSupport::default();
        for &(minus, aligned_length, clip_length) in molecule_support.values() {
            let overlap = aligned_length.min(clip_length);
            support.best_min_overlap = support.best_min_overlap.max(overlap);
            if minus {
                support.minus_reads += 1;
                support.minus_best_min = support.minus_best_min.max(overlap);
            } else {
                support.plus_reads += 1;
                support.plus_best_min = support.plus_best_min.max(overlap);
            }
            support.min_overlap_side1 = match support.min_overlap_side1 {
                0 => aligned_length,
                current => current.min(aligned_length),
            };
            support.min_overlap_side2 = match support.min_overlap_side2 {
                0 => clip_length,
                current => current.min(clip_length),
            };
        }
        let resource_limited = !complete
            || !resource_complete
            || placement_count > MAX_RETAINED_REPEAT_PLACEMENTS
            || family
                .iter()
                .any(|seed| seed.resource_limit_reason.is_some());
        let resolution_is_independent = resolution_molecules
            .iter()
            .any(|molecule_id| !molecule_support.contains_key(molecule_id));
        let state = if resource_limited {
            RepeatAmbiguousState::ResourceLimit
        } else if copies.len() == 1
            && copy_resolved
            && resolution_is_independent
            && anchors_qualified
            && accept_junction(&support)
        {
            RepeatAmbiguousState::Resolved
        } else if copies.len() > 1 && anchors_qualified && accept_junction(&support) {
            RepeatAmbiguousState::SupportedAmbiguous
        } else {
            RepeatAmbiguousState::Candidate
        };
        let feasible_copies: Vec<RepeatCopyPlacement> = copies.into_iter().collect();
        let observed_placements: Vec<RepeatCopyPlacement> = observed.into_iter().collect();
        let breakpoint_end_1 = feasible_copies
            .iter()
            .map(|copy| copy.position_1)
            .max()
            .or_else(|| observed_placements.iter().map(|copy| copy.position_1).max())
            .unwrap_or(first.anchor_position_1);
        let breakpoint_start_1 = feasible_copies
            .iter()
            .map(|copy| copy.position_1)
            .min()
            .or_else(|| observed_placements.iter().map(|copy| copy.position_1).min())
            .unwrap_or(first.anchor_position_1);
        let mut reference_contigs: Vec<String> =
            fasta.iter().map(|record| record.name.clone()).collect();
        reference_contigs.sort();

        RepeatAmbiguousJunctionEvidence {
            state,
            reference_contigs,
            anchor_contig_idx: first.anchor_contig_idx,
            anchor_position_1: first.anchor_position_1,
            anchor_minus: first.anchor_minus,
            feasible_copies,
            observed_placements,
            placement_count,
            placement_family: first.placement_family.clone(),
            placement_sequence: first.placement_sequence.clone(),
            placement_set_complete: complete,
            resource_complete,
            resolution_complete,
            breakpoint_start_1,
            breakpoint_end_1,
            overlap: first.overlap,
            molecule_ids: molecule_support.keys().copied().collect(),
            plus_molecules: support.plus_reads,
            minus_molecules: support.minus_reads,
            clip_lengths: clip_lengths.into_iter().collect(),
            aligned_lengths: aligned_lengths.into_iter().collect(),
            effective_mapqs: effective_mapqs.into_iter().collect(),
            alignment_scores: scores.into_iter().collect(),
            sources: sources.into_iter().collect(),
            pair_geometries: geometry.into_iter().collect(),
            reciprocal_evidence,
            resolution_molecule_ids: resolution_molecules.into_iter().collect(),
            coverage_context: "NOT_USED_FOR_COPY_RESOLUTION",
            rejection_reasons: rejection_reasons.into_iter().collect(),
            exact_breakpoint_supported: state == RepeatAmbiguousState::Resolved,
        }
    }

    pub fn write_tsv(&self, path: &Path, sample_identity: &str) -> crate::error::Result<()> {
        let file = File::create(path)?;
        let mut out = BufWriter::new(file);
        writeln!(out, "sample_id\treference_contigs\tstate\texact_breakpoint_supported\tanchor_contig_idx\tanchor_position_1\tanchor_orientation\tplacement_count\tplacement_set_complete\tfeasible_copy_count\tbreakpoint_interval\toverlap\tdistinct_molecule_count\tplus_molecules\tminus_molecules\tclip_lengths\taligned_lengths\teffective_mapq\talignment_scores\tsources\tpair_geometry\treciprocal_evidence\tresolution_molecule_ids\tcoverage_context\trejection_reasons\tplacement_family\tplacement_sequence\tcandidate_count\trepeat_family_count\toverflow_count\tdiscard_count")?;
        for record in &self.records {
            let mut interval = String::new();
            let _ = write!(
                interval,
                "{}-{}",
                record.breakpoint_start_1, record.breakpoint_end_1
            );
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                sample_identity,
                record.reference_contigs.join(","),
                record.state.as_str(),
                record.exact_breakpoint_supported,
                record.anchor_contig_idx,
                record.anchor_position_1,
                if record.anchor_minus { "-" } else { "+" },
                record.placement_count,
                record.placement_set_complete,
                record.feasible_copies.len(),
                interval,
                record.overlap,
                record.molecule_ids.len(),
                record.plus_molecules,
                record.minus_molecules,
                join_usize(&record.clip_lengths),
                join_usize(&record.aligned_lengths),
                join_u8(&record.effective_mapqs),
                join_i32(&record.alignment_scores),
                join_sources(&record.sources),
                record.pair_geometries.join(","),
                record.reciprocal_evidence,
                join_u64(&record.resolution_molecule_ids),
                record.coverage_context,
                record.rejection_reasons.join(","),
                String::from_utf8_lossy(&record.placement_family),
                String::from_utf8_lossy(&record.placement_sequence),
                self.candidate_count,
                self.repeat_family_count,
                self.overflow_count,
                self.discard_count,
            )?;
        }
        out.flush()?;
        Ok(())
    }
}

fn resolve_reciprocal_seeds(seeds: Vec<RepeatAmbiguousSeed>) -> Vec<RepeatAmbiguousSeed> {
    let mut reverse_index: BTreeMap<RepeatCopyPlacement, Vec<usize>> = BTreeMap::new();
    let mut index_entries = 0usize;
    for (index, seed) in seeds.iter().enumerate() {
        if seed.resource_limit_reason.is_some()
            || !seed.placement_set_complete
            || !seed.resource_complete
            || !seed.unique_anchor_qualified
        {
            continue;
        }
        for placement in &seed.eligible_copies {
            index_entries = index_entries.saturating_add(1);
            if index_entries > MAX_RECIPROCAL_INDEX_ENTRIES {
                return seeds
                    .into_iter()
                    .map(|mut seed| {
                        seed.resource_limit_reason = Some("RESOURCE_LIMIT_RECIPROCAL_INDEX");
                        seed.resource_complete = false;
                        seed.placement_set_complete = false;
                        seed.eligible_copies.clear();
                        seed
                    })
                    .collect();
            }
            reverse_index.entry(*placement).or_default().push(index);
        }
    }

    seeds
        .iter()
        .enumerate()
        .map(|(index, seed)| {
            if seed.resource_limit_reason.is_some()
                || !seed.placement_set_complete
                || !seed.resource_complete
                || !seed.unique_anchor_qualified
                || seed.eligible_copies.len() < 2
            {
                return seed.clone();
            }
            let anchor = seed_anchor(seed);
            let Some(candidate_indexes) = reverse_index.get(&anchor) else {
                return seed.clone();
            };
            if candidate_indexes.len() > MAX_RECIPROCAL_CANDIDATES_PER_SEED {
                let mut limited = seed.clone();
                limited.resource_limit_reason = Some("RESOURCE_LIMIT_RECIPROCAL_CANDIDATES");
                limited.resource_complete = false;
                limited.placement_set_complete = false;
                limited.eligible_copies.clear();
                return limited;
            }

            let evidence = candidate_indexes
                .iter()
                .copied()
                .filter(|other_index| *other_index != index)
                .filter_map(|other_index| {
                    let other = &seeds[other_index];
                    reciprocal_resolution_copy(seed, other).map(|copy| CopyResolutionEvidence {
                        source: CopyResolutionSource::ReciprocalBreakpoint,
                        molecule_id: other.molecule_id,
                        copy,
                    })
                })
                .collect::<Vec<_>>();
            resolve_repeat_copy(seed, &evidence)
        })
        .collect()
}

fn seed_anchor(seed: &RepeatAmbiguousSeed) -> RepeatCopyPlacement {
    RepeatCopyPlacement {
        contig_idx: seed.anchor_contig_idx,
        position_1: seed.anchor_position_1,
        minus: seed.anchor_minus,
    }
}

/// A reciprocal observation is copy-discriminating only when it proves the
/// same two directed endpoints with compatible mosaic geometry. Coordinates
/// alone, low-quality anchors, and re-observation of one molecule are never
/// enough to resolve a repeat copy.
fn reciprocal_resolution_copy(
    seed: &RepeatAmbiguousSeed,
    other: &RepeatAmbiguousSeed,
) -> Option<RepeatCopyPlacement> {
    if seed.molecule_id == other.molecule_id
        || !seed.placement_set_complete
        || !other.placement_set_complete
        || !seed.resource_complete
        || !other.resource_complete
        || seed.resource_limit_reason.is_some()
        || other.resource_limit_reason.is_some()
        || !seed.unique_anchor_qualified
        || !other.unique_anchor_qualified
        || seed.overlap != other.overlap
        || seed.pair_geometry != other.pair_geometry
    {
        return None;
    }
    let other_anchor = seed_anchor(other);
    let copy = seed
        .eligible_copies
        .iter()
        .copied()
        .find(|placement| *placement == other_anchor)?;
    if !other
        .eligible_copies
        .iter()
        .any(|placement| *placement == seed_anchor(seed))
    {
        return None;
    }
    Some(copy)
}

fn join_usize(values: &[usize]) -> String {
    values
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn join_u8(values: &[u8]) -> String {
    values
        .iter()
        .map(u8::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn join_i32(values: &[i32]) -> String {
    values
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn join_u64(values: &[u64]) -> String {
    values
        .iter()
        .map(u64::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn join_sources(values: &[RepeatEvidenceSource]) -> String {
    values
        .iter()
        .map(|source| source.as_str())
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seed(molecule_id: u64, minus: bool, copies: usize) -> RepeatAmbiguousSeed {
        let copies_vec: Vec<RepeatCopyPlacement> = (0..copies)
            .map(|i| RepeatCopyPlacement {
                contig_idx: 0,
                position_1: 200 + i as u64,
                minus: false,
            })
            .collect();
        RepeatAmbiguousSeed {
            source: RepeatEvidenceSource::PrimaryClip,
            anchor_contig_idx: 0,
            anchor_position_1: 100,
            anchor_minus: true,
            observed_placements: copies_vec.clone(),
            eligible_copies: copies_vec,
            placement_count: copies,
            placement_family: b"ACGTACGTACGTAC".to_vec(),
            placement_sequence: b"ACGTACGTACGTAC".to_vec(),
            placement_set_complete: true,
            resource_complete: true,
            resolution_complete: false,
            overlap: 0,
            molecule_id,
            molecule_minus: minus,
            clip_length: 20,
            aligned_length: 20,
            effective_mapq: 42,
            alignment_score: Some(40),
            pair_geometry: "PRIMARY_SOFTCLIP",
            unique_anchor_qualified: true,
            copy_resolved: false,
            resolution_molecule_id: None,
            reciprocal_evidence: false,
            rejection_reason: None,
            resource_limit_reason: None,
        }
    }

    #[test]
    fn supported_multiple_copy_family_is_evidence_only_and_deduplicated() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(
            &fasta,
            vec![seed(7, false, 21), seed(7, false, 21), seed(8, true, 21)],
        );
        let record = &diagnostics.records[0];
        assert_eq!(record.state, RepeatAmbiguousState::SupportedAmbiguous);
        assert!(!record.exact_breakpoint_supported);
        assert_eq!(record.molecule_ids, vec![7, 8]);
        assert_eq!(record.feasible_copies.len(), 21);
    }

    #[test]
    fn two_copy_family_remains_evidence_only() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(
            &fasta,
            vec![seed(7, false, 2), seed(8, true, 2)],
        );
        let record = &diagnostics.records[0];
        assert_eq!(record.state, RepeatAmbiguousState::SupportedAmbiguous);
        assert_eq!(record.feasible_copies.len(), 2);
        assert!(!record.exact_breakpoint_supported);
    }

    #[test]
    fn one_molecule_is_deduplicated_across_primary_and_stage_two_sources() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let primary = seed(7, false, 21);
        let mut stage_two = primary.clone();
        stage_two.source = RepeatEvidenceSource::Stage2Mosaic;
        stage_two.pair_geometry = "STAGE2_MOSAIC";
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, vec![primary, stage_two]);
        let record = &diagnostics.records[0];
        assert_eq!(record.molecule_ids, vec![7]);
        assert_eq!(record.plus_molecules, 1);
        assert_eq!(record.minus_molecules, 0);
    }

    #[test]
    fn singleton_without_independent_copy_resolution_remains_candidate() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(
            &fasta,
            vec![seed(7, false, 1), seed(8, true, 1)],
        );
        assert_eq!(
            diagnostics.records[0].state,
            RepeatAmbiguousState::Candidate
        );
        assert!(!diagnostics.records[0].exact_breakpoint_supported);
    }

    #[test]
    fn independent_copy_resolution_can_mark_a_singleton_resolved() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let mut plus = seed(7, false, 1);
        plus.copy_resolved = true;
        plus.resolution_molecule_id = Some(9);
        let diagnostics =
            RepeatAmbiguousDiagnostics::from_seeds(&fasta, vec![plus, seed(8, true, 1)]);
        assert_eq!(diagnostics.records[0].state, RepeatAmbiguousState::Resolved);
        assert!(diagnostics.records[0].exact_breakpoint_supported);
    }

    #[test]
    fn a_supporting_molecule_cannot_resolve_its_own_copy() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let mut plus = seed(7, false, 1);
        plus.copy_resolved = true;
        plus.resolution_molecule_id = Some(7);
        let diagnostics =
            RepeatAmbiguousDiagnostics::from_seeds(&fasta, vec![plus, seed(8, true, 1)]);
        assert_eq!(
            diagnostics.records[0].state,
            RepeatAmbiguousState::Candidate
        );
        assert!(!diagnostics.records[0].exact_breakpoint_supported);
    }

    #[test]
    fn seed_order_does_not_change_logical_repeat_evidence() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let seeds = vec![seed(8, true, 21), seed(7, false, 21), seed(7, false, 21)];
        let forward = RepeatAmbiguousDiagnostics::from_seeds(&fasta, seeds.clone());
        let reverse =
            RepeatAmbiguousDiagnostics::from_seeds(&fasta, seeds.into_iter().rev().collect());
        assert_eq!(forward, reverse);
    }

    #[test]
    fn unique_mate_geometry_resolves_only_its_supported_copy() {
        let original = seed(7, false, 21);
        let selected = original.eligible_copies[9];
        let resolved = resolve_repeat_copy(
            &original,
            &[CopyResolutionEvidence {
                source: CopyResolutionSource::UniqueMateGeometry,
                molecule_id: 9,
                copy: selected,
            }],
        );
        assert!(resolved.copy_resolved);
        assert_eq!(resolved.eligible_copies, vec![selected]);
        assert_eq!(resolved.placement_count, 21);
    }

    #[test]
    fn conflicting_or_nonindependent_copy_evidence_cannot_resolve() {
        let original = seed(7, false, 21);
        let nonindependent = resolve_repeat_copy(
            &original,
            &[CopyResolutionEvidence {
                source: CopyResolutionSource::UniqueFlankExtension,
                molecule_id: 7,
                copy: original.eligible_copies[0],
            }],
        );
        assert_eq!(nonindependent.eligible_copies.len(), 21);
        let conflicting = resolve_repeat_copy(
            &original,
            &[
                CopyResolutionEvidence {
                    source: CopyResolutionSource::UniqueMateGeometry,
                    molecule_id: 8,
                    copy: original.eligible_copies[0],
                },
                CopyResolutionEvidence {
                    source: CopyResolutionSource::ReciprocalBreakpoint,
                    molecule_id: 9,
                    copy: original.eligible_copies[1],
                },
            ],
        );
        assert_eq!(conflicting.eligible_copies.len(), 21);
    }

    #[test]
    fn distinct_reciprocal_seed_resolves_a_single_repeat_copy() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let forward = seed(7, false, 21);
        let selected = forward.eligible_copies[5];
        let mut reciprocal = seed(9, true, 1);
        reciprocal.anchor_position_1 = selected.position_1;
        reciprocal.anchor_minus = selected.minus;
        reciprocal.observed_placements = vec![RepeatCopyPlacement {
            contig_idx: 0,
            position_1: forward.anchor_position_1,
            minus: forward.anchor_minus,
        }];
        reciprocal.eligible_copies = vec![RepeatCopyPlacement {
            contig_idx: 0,
            position_1: forward.anchor_position_1,
            minus: forward.anchor_minus,
        }];
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(
            &fasta,
            vec![forward, reciprocal, seed(8, true, 21)],
        );
        let resolved = diagnostics
            .records
            .iter()
            .find(|record| record.anchor_position_1 == 100)
            .expect("forward family must be retained");
        assert_eq!(resolved.feasible_copies, vec![selected]);
        assert!(resolved.reciprocal_evidence);
        assert_eq!(resolved.state, RepeatAmbiguousState::Resolved);
        assert!(resolved.exact_breakpoint_supported);
    }

    fn reciprocal_seed(
        forward: &RepeatAmbiguousSeed,
        selected: RepeatCopyPlacement,
        molecule_id: u64,
    ) -> RepeatAmbiguousSeed {
        let mut reciprocal = seed(molecule_id, true, 1);
        reciprocal.anchor_position_1 = selected.position_1;
        reciprocal.anchor_minus = selected.minus;
        reciprocal.observed_placements = vec![seed_anchor(forward)];
        reciprocal.eligible_copies = vec![seed_anchor(forward)];
        reciprocal
    }

    #[test]
    fn reciprocal_resolution_rejects_same_molecule() {
        let forward = seed(7, false, 21);
        let reciprocal = reciprocal_seed(&forward, forward.eligible_copies[5], 7);
        let resolved = resolve_reciprocal_seeds(vec![forward.clone(), reciprocal]);
        assert_eq!(resolved[0].eligible_copies.len(), 21);
        assert!(!resolved[0].copy_resolved);
    }

    #[test]
    fn reciprocal_resolution_rejects_opposite_orientation() {
        let forward = seed(7, false, 21);
        let mut reciprocal = reciprocal_seed(&forward, forward.eligible_copies[5], 9);
        reciprocal.anchor_minus = !reciprocal.anchor_minus;
        let resolved = resolve_reciprocal_seeds(vec![forward.clone(), reciprocal]);
        assert_eq!(resolved[0].eligible_copies.len(), 21);
        assert!(!resolved[0].copy_resolved);
    }

    #[test]
    fn reciprocal_resolution_rejects_unqualified_anchor() {
        let forward = seed(7, false, 21);
        let mut reciprocal = reciprocal_seed(&forward, forward.eligible_copies[5], 9);
        reciprocal.unique_anchor_qualified = false;
        let resolved = resolve_reciprocal_seeds(vec![forward.clone(), reciprocal]);
        assert_eq!(resolved[0].eligible_copies.len(), 21);
        assert!(!resolved[0].copy_resolved);
    }

    #[test]
    fn reciprocal_resolution_rejects_contradictory_geometry() {
        let forward = seed(7, false, 21);
        let mut reciprocal = reciprocal_seed(&forward, forward.eligible_copies[5], 9);
        reciprocal.overlap = 1;
        let resolved = resolve_reciprocal_seeds(vec![forward.clone(), reciprocal]);
        assert_eq!(resolved[0].eligible_copies.len(), 21);
        assert!(!resolved[0].copy_resolved);
    }

    #[test]
    fn reciprocal_resolution_rejects_different_mosaic_geometry() {
        let forward = seed(7, false, 21);
        let mut reciprocal = reciprocal_seed(&forward, forward.eligible_copies[5], 9);
        reciprocal.pair_geometry = "STAGE2_MOSAIC";
        let resolved = resolve_reciprocal_seeds(vec![forward.clone(), reciprocal]);
        assert_eq!(resolved[0].eligible_copies.len(), 21);
        assert!(!resolved[0].copy_resolved);
    }

    #[test]
    fn reciprocal_search_resource_limit_cannot_resolve() {
        let forward = seed(7, false, 21);
        let mut seeds = vec![forward.clone()];
        for molecule_id in 8..=5_008 {
            let mut other = seed(molecule_id, true, 1);
            other.anchor_position_1 = forward.eligible_copies[0].position_1;
            other.anchor_minus = forward.eligible_copies[0].minus;
            other.observed_placements = vec![seed_anchor(&forward)];
            other.eligible_copies = vec![seed_anchor(&forward)];
            seeds.push(other);
        }
        let resolved = resolve_reciprocal_seeds(seeds);
        assert_eq!(
            resolved[0].resource_limit_reason,
            Some("RESOURCE_LIMIT_RECIPROCAL_CANDIDATES")
        );
        assert!(!resolved[0].copy_resolved);
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, resolved);
        assert!(diagnostics
            .records
            .iter()
            .any(|record| record.state == RepeatAmbiguousState::ResourceLimit
                && !record.exact_breakpoint_supported));
    }

    #[test]
    fn diagnostic_sidecar_keeps_state_and_resource_fields() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(
            &fasta,
            vec![seed(7, false, 21), seed(8, true, 21)],
        );
        let path = std::env::temp_dir().join(format!(
            "prokadiff_repeat_ambiguous_{}_{}.tsv",
            std::process::id(),
            diagnostics.records.len()
        ));
        diagnostics.write_tsv(&path, "synthetic-sample").unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        assert!(text.starts_with("sample_id\treference_contigs\tstate\texact_breakpoint_supported"));
        assert!(text.contains("SUPPORTED_AMBIGUOUS"));
        assert!(text.contains("candidate_count"));
        assert!(text.contains("repeat_family_count"));
    }

    #[test]
    fn one_strand_family_remains_candidate() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, vec![seed(7, false, 21)]);
        assert_eq!(
            diagnostics.records[0].state,
            RepeatAmbiguousState::Candidate
        );
    }

    #[test]
    fn incomplete_copy_set_is_resource_limit_not_resolution() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let mut incomplete = seed(7, false, 1);
        incomplete.placement_set_complete = false;
        let diagnostics =
            RepeatAmbiguousDiagnostics::from_seeds(&fasta, vec![incomplete, seed(8, true, 1)]);
        assert_eq!(
            diagnostics.records[0].state,
            RepeatAmbiguousState::ResourceLimit
        );
        assert!(!diagnostics.records[0].exact_breakpoint_supported);
    }

    #[test]
    fn family_scoped_overflow_invalidates_earlier_seeds_in_same_family() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        // Seed 1 is a valid-looking member with 21 eligible copies
        let seed1 = seed(101, false, 21);
        // Seed 2 belongs to the same placement family, but encountered an overflow
        let mut seed2 = seed(102, true, 21);
        seed2.resource_limit_reason = Some("RESOURCE_LIMIT_STAGE2_FAMILY_GEOMETRY");
        seed2.placement_set_complete = false;
        seed2.resource_complete = false;

        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, vec![seed1, seed2]);
        assert_eq!(diagnostics.records.len(), 1);
        let record = &diagnostics.records[0];
        // Whole family must become RESOURCE_LIMIT
        assert_eq!(record.state, RepeatAmbiguousState::ResourceLimit);
        // Feasible copies must be empty (no partial subset may resolve)
        assert!(record.feasible_copies.is_empty());
        assert!(!record.exact_breakpoint_supported);
    }

    #[test]
    fn canonical_molecule_identity_normalizes_flags_and_suffixes() {
        use noodles::sam::alignment::record::Flags;

        // 1. readA/1 with missing mate flag (Flags::empty()) -> mate 1, norm "readA"
        let (norm1, mate1, id1) = canonical_molecule_identity("readA/1", Some(&Flags::empty()));
        assert_eq!(norm1, "readA");
        assert_eq!(mate1, 1);

        // 2. readA with first-mate flag (SEGMENTED | FIRST_SEGMENT) -> mate 1, norm "readA"
        let (norm2, mate2, id2) =
            canonical_molecule_identity("readA", Some(&(Flags::SEGMENTED | Flags::FIRST_SEGMENT)));
        assert_eq!(norm2, "readA");
        assert_eq!(mate2, 1);
        // Cross-path equivalence!
        assert_eq!(id1, id2);

        // 3. readA/2 with missing mate flag -> mate 2, norm "readA"
        let (norm3, mate3, id3) = canonical_molecule_identity("readA/2", Some(&Flags::empty()));
        assert_eq!(norm3, "readA");
        assert_eq!(mate3, 2);

        // 4. readA with second-mate flag (SEGMENTED | LAST_SEGMENT) -> mate 2, norm "readA"
        let (norm4, mate4, id4) =
            canonical_molecule_identity("readA", Some(&(Flags::SEGMENTED | Flags::LAST_SEGMENT)));
        assert_eq!(norm4, "readA");
        assert_eq!(mate4, 2);
        // Cross-path equivalence!
        assert_eq!(id3, id4);

        // Mate 1 and Mate 2 must remain distinguishable!
        assert_ne!(id1, id3);

        // 5. Canonical UNKNOWN mate semantics when flags and suffixes are absent
        let (norm5, mate5, id5) = canonical_molecule_identity("readA", Some(&Flags::empty()));
        assert_eq!(norm5, "readA");
        assert_eq!(mate5, MATE_UNKNOWN);
        assert_ne!(id5, id1);
        assert_ne!(id5, id3);

        let (norm6, mate6, id6) = canonical_molecule_identity("readA", None);
        assert_eq!(norm6, "readA");
        assert_eq!(mate6, MATE_UNKNOWN);
        assert_eq!(id5, id6);
    }

    #[test]
    fn resource_limit_reciprocal_index_triggers_when_exceeding_max_entries() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        // MAX_RECIPROCAL_INDEX_ENTRIES is 200_000.
        // Create seeds whose eligible_copies total > 200_000 (e.g. 50 * 4001 = 200,050).
        let copies_per_seed = 4_001;
        let eligible: Vec<RepeatCopyPlacement> = (0..copies_per_seed)
            .map(|i| RepeatCopyPlacement {
                contig_idx: 0,
                position_1: 100 + i as u64,
                minus: false,
            })
            .collect();
        let mut seeds = Vec::new();
        for mol in 0..50 {
            let mut s = seed(mol, false, copies_per_seed);
            s.anchor_position_1 = 50;
            s.eligible_copies = eligible.clone();
            s.observed_placements = eligible.clone();
            seeds.push(s);
        }
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, seeds);
        assert!(!diagnostics.records.is_empty());
        for rec in &diagnostics.records {
            assert_eq!(rec.state, RepeatAmbiguousState::ResourceLimit);
            assert!(rec
                .rejection_reasons
                .contains(&"RESOURCE_LIMIT_RECIPROCAL_INDEX"));
            assert!(rec.feasible_copies.is_empty());
        }
    }

    #[test]
    fn reordered_reciprocal_inputs_yield_identical_diagnostics() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let forward = seed(7, false, 21);
        let selected = forward.eligible_copies[5];
        let mut reciprocal = seed(9, true, 1);
        reciprocal.anchor_position_1 = selected.position_1;
        reciprocal.anchor_minus = selected.minus;
        reciprocal.observed_placements = vec![RepeatCopyPlacement {
            contig_idx: 0,
            position_1: forward.anchor_position_1,
            minus: forward.anchor_minus,
        }];
        reciprocal.eligible_copies = vec![RepeatCopyPlacement {
            contig_idx: 0,
            position_1: forward.anchor_position_1,
            minus: forward.anchor_minus,
        }];
        let seed_minus = seed(8, true, 21);
        let set_a = vec![forward.clone(), reciprocal.clone(), seed_minus.clone()];
        let set_b = vec![seed_minus, reciprocal, forward];
        let diag_a = RepeatAmbiguousDiagnostics::from_seeds(&fasta, set_a);
        let diag_b = RepeatAmbiguousDiagnostics::from_seeds(&fasta, set_b);
        assert_eq!(diag_a, diag_b);
    }

    #[test]
    fn duplicate_alignment_molecule_deduplication() {
        let fasta = vec![FastaRecord {
            name: "chr".into(),
            seq: vec![b'A'; 1000],
        }];
        let seed1 = seed(42, false, 21);
        let seed2 = seed(42, false, 21);
        let diagnostics = RepeatAmbiguousDiagnostics::from_seeds(&fasta, vec![seed1, seed2]);
        assert_eq!(diagnostics.records.len(), 1);
        assert_eq!(diagnostics.records[0].molecule_ids.len(), 1);
        assert_eq!(diagnostics.records[0].molecule_ids[0], 42);
    }
}
