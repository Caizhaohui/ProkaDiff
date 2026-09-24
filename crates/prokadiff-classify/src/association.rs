use prokadiff_gd::GdKind;
use prokadiff_offtarget::{BulgeType, MutationOffTargetLink, OffTargetSite, Strand};
use std::collections::BTreeMap;

use crate::differential::{DifferentialEvent, EventId};
use crate::geometry::{
    event_geometry, interval_distance, junction_side_distance, JunctionSideTag, MutationGeometry,
};

pub const DEFAULT_NEAR_DISTANCE: u64 = 50;
pub const DEFAULT_MAX_MISMATCHES: u32 = 4;

/// Status of the candidate off-target search.
///
/// Disambiguates search not requested/performed from search performed finding zero candidates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CandidateSearchStatus {
    NotPerformed,
    PerformedNoCandidates,
    PerformedWithCandidates,
}

impl CandidateSearchStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NotPerformed => "NOT_PERFORMED",
            Self::PerformedNoCandidates => "PERFORMED_NO_CANDIDATES",
            Self::PerformedWithCandidates => "PERFORMED_WITH_CANDIDATES",
        }
    }
}

/// Unified candidate off-target search outcome.
#[derive(Clone, Debug, PartialEq)]
pub struct CandidateSearchOutcome {
    pub status: CandidateSearchStatus,
    pub sites: Vec<OffTargetSite>,
}

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum AssociationProjectionError {
    #[error("association references unknown differential event {0}")]
    UnknownEventId(EventId),
}

impl Default for CandidateSearchOutcome {
    fn default() -> Self {
        Self {
            status: CandidateSearchStatus::NotPerformed,
            sites: Vec::new(),
        }
    }
}

/// Spatial and sequence association between a differential mutation event and a candidate off-target site.
///
/// Preserves evidence tier (backend, bulge type/size) and observed genomic PAM.
/// Spatial association represents proximity within the specified window, not validated causal cleavage.
#[derive(Clone, Debug, PartialEq)]
pub struct MutationSiteAssociation {
    pub event_id: EventId,
    pub mutation_type: String,
    pub mutation_position: u64,
    pub junction_side: Option<JunctionSideTag>,
    pub site_id: String,
    pub site_seq_id: String,
    pub site_start: u64,
    pub site_end: u64,
    pub site_strand: Strand,
    pub distance_to_site: u64,
    pub mismatches: u32,
    pub pam: String,
    pub cfd_score: Option<f64>,
    pub hsu_score: Option<f64>,
    pub search_backend: String,
    pub bulge_type: BulgeType,
    pub bulge_size: u32,
    pub target_seq: String,
    pub association_window: u64,
}

/// Sort candidate off-target sites deterministically according to total ordering and renumber site IDs.
pub fn canonical_sort_sites(sites: &mut [OffTargetSite]) {
    sites.sort_by(|a, b| {
        a.seq_id
            .cmp(&b.seq_id)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
            .then_with(|| a.strand.cmp(&b.strand))
            .then_with(|| a.mismatches.cmp(&b.mismatches))
            .then_with(|| a.target_seq.cmp(&b.target_seq))
            .then_with(|| a.pam.cmp(&b.pam))
    });
    for (idx, site) in sites.iter_mut().enumerate() {
        site.site_id = format!("SITE_{:06}", idx + 1);
    }
}

/// Sort mutation site associations deterministically.
///
/// Unvalidated scores (cfd_score, hsu_score) remain metadata only and do not affect ordering.
pub fn canonical_sort_associations(assocs: &mut [MutationSiteAssociation]) {
    assocs.sort_by(|a, b| {
        a.event_id
            .cmp(&b.event_id)
            .then_with(|| a.junction_side.cmp(&b.junction_side))
            .then_with(|| a.distance_to_site.cmp(&b.distance_to_site))
            .then_with(|| a.mismatches.cmp(&b.mismatches))
            .then_with(|| a.site_id.cmp(&b.site_id))
            .then_with(|| a.mutation_position.cmp(&b.mutation_position))
    });
}

/// Deterministic, non-destructive read-time projection of primary association for an event.
///
/// Tie-breaking order:
/// 1. Minimum distance_to_site
/// 2. Minimum mismatches
/// 3. Junction side (Side1 < Side2)
/// 4. Lexicographical site_id
///
/// Note: cfd_score and hsu_score do NOT affect ordering; scoring remains metadata only
/// while validation gates are disabled.
pub fn primary_association(assocs: &[MutationSiteAssociation]) -> Option<&MutationSiteAssociation> {
    assocs.iter().min_by(|a, b| {
        a.distance_to_site
            .cmp(&b.distance_to_site)
            .then_with(|| a.mismatches.cmp(&b.mismatches))
            .then_with(|| a.junction_side.cmp(&b.junction_side))
            .then_with(|| a.site_id.cmp(&b.site_id))
    })
}

fn build_association(
    event_id: &EventId,
    mutation_type: &str,
    mutation_position: u64,
    junction_side: Option<JunctionSideTag>,
    site: &OffTargetSite,
    distance_to_site: u64,
    association_window: u64,
) -> MutationSiteAssociation {
    MutationSiteAssociation {
        event_id: event_id.clone(),
        mutation_type: mutation_type.to_string(),
        mutation_position,
        junction_side,
        site_id: site.site_id.clone(),
        site_seq_id: site.seq_id.clone(),
        site_start: site.start,
        site_end: site.end,
        site_strand: site.strand,
        distance_to_site,
        mismatches: site.mismatches,
        pam: site.pam.clone(),
        cfd_score: site.cfd_score,
        hsu_score: site.hsu_score,
        search_backend: site.search_backend.clone(),
        bulge_type: site.bulge_type,
        bulge_size: site.bulge_size,
        target_seq: site.target_seq.clone(),
        association_window,
    }
}

/// Linear reference implementation for associating differential events to off-target sites.
///
/// Used for correctness verification and gold-standard oracle checks against the indexed implementation.
pub fn linear_associate_events_to_sites(
    events: &[DifferentialEvent],
    sites: &[OffTargetSite],
    window_bp: u64,
) -> Vec<MutationSiteAssociation> {
    let mut assocs = Vec::new();
    for event in events {
        let entry = &event.representative;
        // Amendment 2: MC, RA, and UN evidence-only records must never produce MutationSiteAssociation
        if matches!(entry.kind, GdKind::Mc | GdKind::Ra | GdKind::Un) {
            continue;
        }
        let Some(geom) = event_geometry(entry) else {
            continue;
        };

        match geom {
            MutationGeometry::Point { seq_id, position } => {
                for site in sites {
                    if site.seq_id != seq_id {
                        continue;
                    }
                    let dist = interval_distance(position, position, site.start, site.end);
                    if dist <= window_bp {
                        assocs.push(build_association(
                            &event.event_id,
                            entry.kind.as_str(),
                            position,
                            None,
                            site,
                            dist,
                            window_bp,
                        ));
                    }
                }
            }
            MutationGeometry::Span { seq_id, start, end } => {
                let pos = entry.position().unwrap_or(start);
                for site in sites {
                    if site.seq_id != seq_id {
                        continue;
                    }
                    let dist = interval_distance(start, end, site.start, site.end);
                    if dist <= window_bp {
                        assocs.push(build_association(
                            &event.event_id,
                            entry.kind.as_str(),
                            pos,
                            None,
                            site,
                            dist,
                            window_bp,
                        ));
                    }
                }
            }
            MutationGeometry::Junction { side1, side2 } => {
                // Amendment 3: Preserve ALL qualifying breakpoint-site associations.
                for site in sites {
                    if let Some(dist) = junction_side_distance(&side1, site) {
                        if dist <= window_bp {
                            assocs.push(build_association(
                                &event.event_id,
                                entry.kind.as_str(),
                                side1.position,
                                Some(JunctionSideTag::Side1),
                                site,
                                dist,
                                window_bp,
                            ));
                        }
                    }
                }
                for site in sites {
                    if let Some(dist) = junction_side_distance(&side2, site) {
                        if dist <= window_bp {
                            assocs.push(build_association(
                                &event.event_id,
                                entry.kind.as_str(),
                                side2.position,
                                Some(JunctionSideTag::Side2),
                                site,
                                dist,
                                window_bp,
                            ));
                        }
                    }
                }
            }
        }
    }
    canonical_sort_associations(&mut assocs);
    assocs
}

struct ContigIndex<'a> {
    sites: Vec<&'a OffTargetSite>,
    max_site_len: u64,
}

impl<'a> ContigIndex<'a> {
    fn query_sites(&self, q_start: u64, q_end: u64, window_bp: u64) -> &[&'a OffTargetSite] {
        let min_start = q_start.saturating_sub(window_bp.saturating_add(self.max_site_len));
        let max_start = q_end.saturating_add(window_bp);
        let low_idx = self.sites.partition_point(|s| s.start < min_start);
        let high_idx = self.sites.partition_point(|s| s.start <= max_start);
        &self.sites[low_idx..high_idx]
    }
}

/// Associating differential events to off-target sites using coordinate-sorted indexing and window restriction.
///
/// Provably equivalent to linear scan while optimizing lookup complexity.
pub fn associate_events_to_sites(
    events: &[DifferentialEvent],
    sites: &[OffTargetSite],
    window_bp: u64,
) -> Vec<MutationSiteAssociation> {
    if events.is_empty() || sites.is_empty() {
        return Vec::new();
    }

    // Build per-contig coordinate-sorted index
    let mut contig_map: BTreeMap<&str, Vec<&OffTargetSite>> = BTreeMap::new();
    for site in sites {
        contig_map.entry(&site.seq_id).or_default().push(site);
    }

    let mut indexed_contigs: BTreeMap<&str, ContigIndex> = BTreeMap::new();
    for (seq_id, mut c_sites) in contig_map {
        c_sites.sort_by_key(|s| s.start);
        let max_site_len = c_sites
            .iter()
            .map(|s| s.end.saturating_sub(s.start).saturating_add(1))
            .max()
            .unwrap_or(100)
            .max(100);
        indexed_contigs.insert(
            seq_id,
            ContigIndex {
                sites: c_sites,
                max_site_len,
            },
        );
    }

    let mut assocs = Vec::new();

    for event in events {
        let entry = &event.representative;
        // Amendment 2: MC, RA, and UN evidence-only records must never produce MutationSiteAssociation
        if matches!(entry.kind, GdKind::Mc | GdKind::Ra | GdKind::Un) {
            continue;
        }
        let Some(geom) = event_geometry(entry) else {
            continue;
        };

        match geom {
            MutationGeometry::Point { seq_id, position } => {
                if let Some(c_index) = indexed_contigs.get(seq_id.as_str()) {
                    for site in c_index.query_sites(position, position, window_bp) {
                        let dist = interval_distance(position, position, site.start, site.end);
                        if dist <= window_bp {
                            assocs.push(build_association(
                                &event.event_id,
                                entry.kind.as_str(),
                                position,
                                None,
                                site,
                                dist,
                                window_bp,
                            ));
                        }
                    }
                }
            }
            MutationGeometry::Span { seq_id, start, end } => {
                let pos = entry.position().unwrap_or(start);
                if let Some(c_index) = indexed_contigs.get(seq_id.as_str()) {
                    for site in c_index.query_sites(start, end, window_bp) {
                        let dist = interval_distance(start, end, site.start, site.end);
                        if dist <= window_bp {
                            assocs.push(build_association(
                                &event.event_id,
                                entry.kind.as_str(),
                                pos,
                                None,
                                site,
                                dist,
                                window_bp,
                            ));
                        }
                    }
                }
            }
            MutationGeometry::Junction { side1, side2 } => {
                // Amendment 3: Preserve ALL qualifying breakpoint-site associations.
                if let Some(c_index) = indexed_contigs.get(side1.seq_id.as_str()) {
                    for site in c_index.query_sites(side1.position, side1.position, window_bp) {
                        if let Some(dist) = junction_side_distance(&side1, site) {
                            if dist <= window_bp {
                                assocs.push(build_association(
                                    &event.event_id,
                                    entry.kind.as_str(),
                                    side1.position,
                                    Some(JunctionSideTag::Side1),
                                    site,
                                    dist,
                                    window_bp,
                                ));
                            }
                        }
                    }
                }
                if let Some(c_index) = indexed_contigs.get(side2.seq_id.as_str()) {
                    for site in c_index.query_sites(side2.position, side2.position, window_bp) {
                        if let Some(dist) = junction_side_distance(&side2, site) {
                            if dist <= window_bp {
                                assocs.push(build_association(
                                    &event.event_id,
                                    entry.kind.as_str(),
                                    side2.position,
                                    Some(JunctionSideTag::Side2),
                                    site,
                                    dist,
                                    window_bp,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    canonical_sort_associations(&mut assocs);
    assocs
}

/// Project internal EventId-based `MutationSiteAssociation` records to public `MutationOffTargetLink` TSV format.
///
/// Preserves legacy `mut_<id>` identifier for M3 public compatibility while mapping from canonical EventId.
pub fn project_associations_to_links(
    associations: &[MutationSiteAssociation],
    events: &[DifferentialEvent],
) -> Result<Vec<MutationOffTargetLink>, AssociationProjectionError> {
    associations
        .iter()
        .map(|assoc| {
            let mut_id = events
                .iter()
                .find(|e| e.event_id == assoc.event_id)
                .map(|e| format!("mut_{}", e.representative.id))
                .ok_or_else(|| {
                    AssociationProjectionError::UnknownEventId(assoc.event_id.clone())
                })?;
            Ok(MutationOffTargetLink {
                mutation_id: mut_id,
                site_id: assoc.site_id.clone(),
                mutation_type: assoc.mutation_type.clone(),
                mutation_position: assoc.mutation_position,
                site_start: assoc.site_start,
                site_end: assoc.site_end,
                distance_to_site: assoc.distance_to_site,
                mismatches: assoc.mismatches,
                pam: assoc.pam.clone(),
                cfd_score: assoc.cfd_score,
                association_window: assoc.association_window,
            })
        })
        .collect()
}
