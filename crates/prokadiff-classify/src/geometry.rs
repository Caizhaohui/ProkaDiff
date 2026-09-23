use prokadiff_gd::{GdEntry, GdKind};
use prokadiff_offtarget::{OffTargetSite, Strand};
use std::str::FromStr;

use crate::RefContig;

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum GeometryError {
    #[error("entry {kind} is evidence-only or lacks genomic coordinates", kind = .0.as_str())]
    NoGeometry(GdKind),
    #[error("invalid 1-based coordinate: position cannot be 0")]
    ZeroCoordinate,
    #[error("invalid coordinate span: start {start} > end {end}")]
    InvertedSpan { start: u64, end: u64 },
    #[error("coordinate overflow: start {start} + size {size} exceeds u64::MAX")]
    CoordinateOverflow { start: u64, size: u64 },
    #[error("coordinate {end} exceeds contig '{seq_id}' boundary of length {contig_len}")]
    ContigBoundaryExceeded {
        seq_id: String,
        end: u64,
        contig_len: u64,
    },
    #[error("unknown contig '{0}' not found in reference sequences")]
    UnknownContig(String),
    #[error("missing required field in entry: {0}")]
    MissingField(String),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum JunctionSideTag {
    Side1,
    Side2,
}

impl JunctionSideTag {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Side1 => "side_1",
            Self::Side2 => "side_2",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct JunctionSide {
    pub side: JunctionSideTag,
    pub seq_id: String,
    pub position: u64,
    pub strand: Option<Strand>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum MutationGeometry {
    Point {
        seq_id: String,
        position: u64,
    },
    Span {
        seq_id: String,
        start: u64,
        end: u64,
    },
    Junction {
        side1: JunctionSide,
        side2: JunctionSide,
    },
}

impl MutationGeometry {
    /// Return the list of closed 1-based intervals `(seq_id, start, end)` representing this geometry.
    /// Exactly preserves entry_intervals() coordinate semantics.
    pub fn intervals(&self) -> Vec<(String, u64, u64)> {
        match self {
            Self::Point { seq_id, position } => vec![(seq_id.clone(), *position, *position)],
            Self::Span { seq_id, start, end } => vec![(seq_id.clone(), *start, *end)],
            Self::Junction { side1, side2 } => vec![
                (side1.seq_id.clone(), side1.position, side1.position),
                (side2.seq_id.clone(), side2.position, side2.position),
            ],
        }
    }

    /// Primary contig / sequence ID for single-locus mutations (or side 1 for junctions).
    pub fn primary_seq_id(&self) -> &str {
        match self {
            Self::Point { seq_id, .. } | Self::Span { seq_id, .. } => seq_id.as_str(),
            Self::Junction { side1, .. } => side1.seq_id.as_str(),
        }
    }

    /// Primary 1-based start coordinate.
    pub fn primary_start(&self) -> u64 {
        match self {
            Self::Point { position, .. } => *position,
            Self::Span { start, .. } => *start,
            Self::Junction { side1, .. } => side1.position,
        }
    }

    /// Primary 1-based end coordinate.
    pub fn primary_end(&self) -> u64 {
        match self {
            Self::Point { position, .. } => *position,
            Self::Span { end, .. } => *end,
            Self::Junction { side1, .. } => side1.position,
        }
    }
}

/// Parse geometry from a Genome Diff entry.
///
/// Exactly preserves existing entry_intervals() semantics while strictly validating coordinates.
pub fn event_geometry(e: &GdEntry) -> Option<MutationGeometry> {
    match e.kind {
        GdKind::Ra => None,
        GdKind::Un | GdKind::Mc => {
            let seq = e.fields.first()?.clone();
            let start = e.fields.get(1)?.parse::<u64>().ok()?;
            let end = e.fields.get(2)?.parse::<u64>().ok()?;
            if start == 0 || end < start {
                return None;
            }
            Some(MutationGeometry::Span {
                seq_id: seq,
                start,
                end,
            })
        }
        GdKind::Jc => {
            let s1 = e.fields.first()?.clone();
            let p1 = e.fields.get(1)?.parse::<u64>().ok()?;
            let st1 = e.fields.get(2).and_then(|s| Strand::from_str(s).ok());

            let s2 = e.fields.get(3)?.clone();
            let p2 = e.fields.get(4)?.parse::<u64>().ok()?;
            let st2 = e.fields.get(5).and_then(|s| Strand::from_str(s).ok());

            if p1 == 0 || p2 == 0 {
                return None;
            }

            Some(MutationGeometry::Junction {
                side1: JunctionSide {
                    side: JunctionSideTag::Side1,
                    seq_id: s1,
                    position: p1,
                    strand: st1,
                },
                side2: JunctionSide {
                    side: JunctionSideTag::Side2,
                    seq_id: s2,
                    position: p2,
                    strand: st2,
                },
            })
        }
        GdKind::Del | GdKind::Sub | GdKind::Inv | GdKind::Amp | GdKind::Con => {
            let seq = e.fields.first()?.clone();
            let start = e.fields.get(1)?.parse::<u64>().ok()?;
            if start == 0 {
                return None;
            }
            let size = e
                .fields
                .get(2)
                .and_then(|x| x.parse::<u64>().ok())
                .unwrap_or(1)
                .max(1);
            let end = start.checked_add(size - 1)?;
            Some(MutationGeometry::Span {
                seq_id: seq,
                start,
                end,
            })
        }
        _ => {
            let seq = e.seq_id()?.to_string();
            let pos = e.position()?;
            if pos == 0 {
                return None;
            }
            Some(MutationGeometry::Point {
                seq_id: seq,
                position: pos,
            })
        }
    }
}

/// Parse geometry from a Genome Diff entry with validation against reference contig boundaries.
pub fn event_geometry_checked(
    e: &GdEntry,
    refs: &[RefContig],
) -> Result<MutationGeometry, GeometryError> {
    let geom = match e.kind {
        GdKind::Ra => return Err(GeometryError::NoGeometry(e.kind)),
        GdKind::Un | GdKind::Mc => {
            let seq = e
                .fields
                .first()
                .cloned()
                .ok_or_else(|| GeometryError::MissingField("seq_id".into()))?;
            let start = e
                .fields
                .get(1)
                .and_then(|x| x.parse::<u64>().ok())
                .ok_or_else(|| GeometryError::MissingField("start".into()))?;
            let end = e
                .fields
                .get(2)
                .and_then(|x| x.parse::<u64>().ok())
                .ok_or_else(|| GeometryError::MissingField("end".into()))?;
            if start == 0 {
                return Err(GeometryError::ZeroCoordinate);
            }
            if end < start {
                return Err(GeometryError::InvertedSpan { start, end });
            }
            MutationGeometry::Span {
                seq_id: seq,
                start,
                end,
            }
        }
        GdKind::Jc => {
            let s1 = e
                .fields
                .first()
                .cloned()
                .ok_or_else(|| GeometryError::MissingField("side1_seq_id".into()))?;
            let p1 = e
                .fields
                .get(1)
                .and_then(|x| x.parse::<u64>().ok())
                .ok_or_else(|| GeometryError::MissingField("side1_position".into()))?;
            let st1 = e.fields.get(2).and_then(|s| Strand::from_str(s).ok());

            let s2 = e
                .fields
                .get(3)
                .cloned()
                .ok_or_else(|| GeometryError::MissingField("side2_seq_id".into()))?;
            let p2 = e
                .fields
                .get(4)
                .and_then(|x| x.parse::<u64>().ok())
                .ok_or_else(|| GeometryError::MissingField("side2_position".into()))?;
            let st2 = e.fields.get(5).and_then(|s| Strand::from_str(s).ok());

            if p1 == 0 || p2 == 0 {
                return Err(GeometryError::ZeroCoordinate);
            }

            MutationGeometry::Junction {
                side1: JunctionSide {
                    side: JunctionSideTag::Side1,
                    seq_id: s1,
                    position: p1,
                    strand: st1,
                },
                side2: JunctionSide {
                    side: JunctionSideTag::Side2,
                    seq_id: s2,
                    position: p2,
                    strand: st2,
                },
            }
        }
        GdKind::Del | GdKind::Sub | GdKind::Inv | GdKind::Amp | GdKind::Con => {
            let seq = e
                .fields
                .first()
                .cloned()
                .ok_or_else(|| GeometryError::MissingField("seq_id".into()))?;
            let start = e
                .fields
                .get(1)
                .and_then(|x| x.parse::<u64>().ok())
                .ok_or_else(|| GeometryError::MissingField("position".into()))?;
            if start == 0 {
                return Err(GeometryError::ZeroCoordinate);
            }
            let size = e
                .fields
                .get(2)
                .and_then(|x| x.parse::<u64>().ok())
                .unwrap_or(1)
                .max(1);
            let end = start
                .checked_add(size - 1)
                .ok_or(GeometryError::CoordinateOverflow { start, size })?;
            MutationGeometry::Span {
                seq_id: seq,
                start,
                end,
            }
        }
        _ => {
            let seq = e
                .seq_id()
                .map(str::to_string)
                .ok_or_else(|| GeometryError::MissingField("seq_id".into()))?;
            let pos = e
                .position()
                .ok_or_else(|| GeometryError::MissingField("position".into()))?;
            if pos == 0 {
                return Err(GeometryError::ZeroCoordinate);
            }
            MutationGeometry::Point {
                seq_id: seq,
                position: pos,
            }
        }
    };

    // Check contig boundary limits
    for (sid, start, end) in geom.intervals() {
        if let Some(contig) = refs.iter().find(|r| r.name == sid) {
            let len = contig.seq.len() as u64;
            if start > len || end > len {
                return Err(GeometryError::ContigBoundaryExceeded {
                    seq_id: sid,
                    end: end.max(start),
                    contig_len: len,
                });
            }
        } else {
            return Err(GeometryError::UnknownContig(sid));
        }
    }

    Ok(geom)
}

/// Helper returning all geometries parsed for an entry (single element or empty for RA).
pub fn event_geometries(e: &GdEntry) -> Vec<MutationGeometry> {
    match event_geometry(e) {
        Some(g) => vec![g],
        None => Vec::new(),
    }
}

/// Calculate distance between two 1-based closed genomic intervals `[a_start, a_end]` and `[b_start, b_end]`.
///
/// - Overlapping intervals return 0 bp.
/// - Adjacent intervals (e.g. 100..100 and 101..101) return 1 bp.
/// - Returns 0 if coordinates overlap.
pub fn interval_distance(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> u64 {
    if a_start <= b_end && b_start <= a_end {
        0
    } else if a_start > b_end {
        a_start - b_end
    } else {
        b_start - a_end
    }
}

/// Calculate distance between a 1-based point coordinate and a 1-based closed target interval.
pub fn point_distance(pos: u64, b_start: u64, b_end: u64) -> u64 {
    interval_distance(pos, pos, b_start, b_end)
}

/// Calculate the minimum distance between a MutationGeometry and an OffTargetSite.
///
/// Returns None if the mutation and site are on different contigs/chromosomes.
pub fn geometry_distance(geom: &MutationGeometry, site: &OffTargetSite) -> Option<u64> {
    match geom {
        MutationGeometry::Point { seq_id, position } => {
            if seq_id != &site.seq_id {
                None
            } else {
                Some(point_distance(*position, site.start, site.end))
            }
        }
        MutationGeometry::Span { seq_id, start, end } => {
            if seq_id != &site.seq_id {
                None
            } else {
                Some(interval_distance(*start, *end, site.start, site.end))
            }
        }
        MutationGeometry::Junction { side1, side2 } => {
            let d1 = if side1.seq_id == site.seq_id {
                Some(point_distance(side1.position, site.start, site.end))
            } else {
                None
            };
            let d2 = if side2.seq_id == site.seq_id {
                Some(point_distance(side2.position, site.start, site.end))
            } else {
                None
            };
            match (d1, d2) {
                (Some(a), Some(b)) => Some(a.min(b)),
                (Some(a), None) => Some(a),
                (None, Some(b)) => Some(b),
                (None, None) => None,
            }
        }
    }
}

/// Calculate distance between a specific JunctionSide breakpoint and an OffTargetSite.
pub fn junction_side_distance(side: &JunctionSide, site: &OffTargetSite) -> Option<u64> {
    if side.seq_id != site.seq_id {
        None
    } else {
        Some(point_distance(side.position, site.start, site.end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mock_refs() -> Vec<RefContig> {
        vec![
            RefContig {
                name: "chr1".into(),
                seq: vec![b'A'; 1000],
            },
            RefContig {
                name: "chr2".into(),
                seq: vec![b'A'; 2000],
            },
        ]
    }

    #[test]
    fn test_point_geometry_and_intervals() {
        let snp = GdEntry::snp(1, "chr1", 100, "T");
        let geom = event_geometry(&snp).expect("valid point geometry");
        assert_eq!(
            geom,
            MutationGeometry::Point {
                seq_id: "chr1".into(),
                position: 100,
            }
        );
        assert_eq!(geom.intervals(), vec![("chr1".to_string(), 100, 100)]);
    }

    #[test]
    fn test_size_1_span() {
        let del1 = GdEntry::del(2, "chr1", 50, 1);
        let geom = event_geometry(&del1).expect("valid size-1 deletion");
        assert_eq!(
            geom,
            MutationGeometry::Span {
                seq_id: "chr1".into(),
                start: 50,
                end: 50,
            }
        );
        assert_eq!(geom.intervals(), vec![("chr1".to_string(), 50, 50)]);
    }

    #[test]
    fn test_normal_del_sub_inv_amp_con_span() {
        // DEL: start 100, size 50 -> 100..=149
        let del = GdEntry::del(3, "chr1", 100, 50);
        assert_eq!(
            event_geometry(&del).unwrap().intervals(),
            vec![("chr1".to_string(), 100, 149)]
        );

        // SUB: start 200, size 30, alt "..." -> 200..=229
        let sub = GdEntry::sub(4, "chr1", 200, 30, "ACGT");
        assert_eq!(
            event_geometry(&sub).unwrap().intervals(),
            vec![("chr1".to_string(), 200, 229)]
        );

        // INV: start 300, size 100 -> 300..=399
        let inv = GdEntry::inv(5, "chr1", 300, 100);
        assert_eq!(
            event_geometry(&inv).unwrap().intervals(),
            vec![("chr1".to_string(), 300, 399)]
        );

        // AMP: start 400, size 25, copy_num 2 -> 400..=424
        let amp = GdEntry::amp(6, "chr1", 400, 25, 2);
        assert_eq!(
            event_geometry(&amp).unwrap().intervals(),
            vec![("chr1".to_string(), 400, 424)]
        );

        // CON: start 500, size 40, region "chr1:1-40" -> 500..=539
        let con = GdEntry::con(7, "chr1", 500, 40, "chr1:1-40");
        assert_eq!(
            event_geometry(&con).unwrap().intervals(),
            vec![("chr1".to_string(), 500, 539)]
        );
    }

    #[test]
    fn test_span_ending_exactly_at_contig_boundary() {
        let refs = mock_refs(); // chr1 len = 1000

        // Span 901..=1000 (size 100) exactly hits boundary
        let del = GdEntry::del(10, "chr1", 901, 100);
        let checked = event_geometry_checked(&del, &refs);
        assert!(checked.is_ok());
        assert_eq!(
            checked.unwrap(),
            MutationGeometry::Span {
                seq_id: "chr1".into(),
                start: 901,
                end: 1000,
            }
        );

        // Span 901..=1001 (size 101) exceeds boundary
        let del_overflow = GdEntry::del(11, "chr1", 901, 101);
        let err = event_geometry_checked(&del_overflow, &refs);
        assert_eq!(
            err,
            Err(GeometryError::ContigBoundaryExceeded {
                seq_id: "chr1".into(),
                end: 1001,
                contig_len: 1000,
            })
        );
    }

    #[test]
    fn test_checked_overflow_and_invalid_coordinates() {
        let refs = mock_refs();

        // Position 0 is rejected
        let zero_del = GdEntry::del(20, "chr1", 0, 10);
        assert_eq!(event_geometry(&zero_del), None);
        assert_eq!(
            event_geometry_checked(&zero_del, &refs),
            Err(GeometryError::ZeroCoordinate)
        );

        // u64 overflow in start + size - 1
        let mut overflow_del = GdEntry::del(21, "chr1", u64::MAX - 5, 10);
        overflow_del.fields[1] = (u64::MAX - 5).to_string();
        overflow_del.fields[2] = "10".to_string();
        assert_eq!(event_geometry(&overflow_del), None);
        assert_eq!(
            event_geometry_checked(&overflow_del, &refs),
            Err(GeometryError::CoordinateOverflow {
                start: u64::MAX - 5,
                size: 10,
            })
        );
    }

    #[test]
    fn test_jc_geometry() {
        let jc = GdEntry::jc(30, "chr1", 100, "+", "chr2", 500, "-", 0);
        let geom = event_geometry(&jc).expect("valid junction");
        assert_eq!(
            geom,
            MutationGeometry::Junction {
                side1: JunctionSide {
                    side: JunctionSideTag::Side1,
                    seq_id: "chr1".into(),
                    position: 100,
                    strand: Some(Strand::Plus),
                },
                side2: JunctionSide {
                    side: JunctionSideTag::Side2,
                    seq_id: "chr2".into(),
                    position: 500,
                    strand: Some(Strand::Minus),
                },
            }
        );
        assert_eq!(
            geom.intervals(),
            vec![
                ("chr1".to_string(), 100, 100),
                ("chr2".to_string(), 500, 500),
            ]
        );
    }

    #[test]
    fn test_evidence_only_ra_has_no_geometry() {
        let ra = GdEntry::ra(40, "chr1", 100, 0, "A", "T");
        assert_eq!(event_geometry(&ra), None);
        assert!(event_geometries(&ra).is_empty());
    }

    #[test]
    fn test_checked_unknown_contig_returns_typed_error() {
        let refs = mock_refs();
        let unknown_snp = GdEntry::snp(50, "unknown_chr", 100, "C");
        assert_eq!(
            event_geometry_checked(&unknown_snp, &refs),
            Err(GeometryError::UnknownContig("unknown_chr".into()))
        );

        let unknown_jc = GdEntry::jc(51, "chr1", 100, "+", "unknown_chr2", 50, "-", 0);
        assert_eq!(
            event_geometry_checked(&unknown_jc, &refs),
            Err(GeometryError::UnknownContig("unknown_chr2".into()))
        );
    }

    #[test]
    fn test_checked_known_contig_valid_coordinate() {
        let refs = mock_refs();
        let valid_snp = GdEntry::snp(52, "chr1", 100, "C");
        assert_eq!(
            event_geometry_checked(&valid_snp, &refs),
            Ok(MutationGeometry::Point {
                seq_id: "chr1".into(),
                position: 100,
            })
        );

        let valid_del = GdEntry::del(53, "chr1", 100, 50);
        assert_eq!(
            event_geometry_checked(&valid_del, &refs),
            Ok(MutationGeometry::Span {
                seq_id: "chr1".into(),
                start: 100,
                end: 149,
            })
        );
    }

    #[test]
    fn test_checked_known_contig_out_of_range_coordinate() {
        let refs = mock_refs(); // chr1 len = 1000
        let oob_snp = GdEntry::snp(54, "chr1", 1001, "C");
        assert_eq!(
            event_geometry_checked(&oob_snp, &refs),
            Err(GeometryError::ContigBoundaryExceeded {
                seq_id: "chr1".into(),
                end: 1001,
                contig_len: 1000,
            })
        );

        let oob_del = GdEntry::del(55, "chr1", 950, 60); // 950..=1009 > 1000
        assert_eq!(
            event_geometry_checked(&oob_del, &refs),
            Err(GeometryError::ContigBoundaryExceeded {
                seq_id: "chr1".into(),
                end: 1009,
                contig_len: 1000,
            })
        );
    }
}
