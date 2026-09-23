//! Canonical differential event representation, deterministic DV1 event IDs,
//! and evidence lineage for ProkaDiff (M1).

use prokadiff_gd::{right_align_del, right_align_ins, GdEntry, GdKind, GenomeDiff};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashSet};
use std::fmt;
use thiserror::Error;

use crate::{is_product_mutation, RefContig};

mod matching;

use matching::{minimum_cost_maximum_matching, CoordinateCost, TolerantCandidate};

const STRUCTURAL_TOLERANCE_BP: u64 = 5;
const DNA_IUPAC_ALPHABET: &[u8] = b"ACGTRYSWKMBDHVN";

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum DifferentialError {
    #[error("contig '{0}' not found in reference skeleton")]
    ContigNotFound(String),

    #[error("invalid sequence for contig '{seq_id}': {msg}")]
    InvalidContigSequence { seq_id: String, msg: String },

    #[error("invalid 1-based position {position} for contig '{seq_id}' (length {contig_len})")]
    InvalidPosition {
        seq_id: String,
        position: u64,
        contig_len: usize,
    },

    #[error("invalid size {size} for contig '{seq_id}' at position {position} (contig length {contig_len})")]
    InvalidSize {
        seq_id: String,
        position: u64,
        size: u64,
        contig_len: usize,
    },

    #[error("invalid allele '{allele}' for contig '{seq_id}'")]
    InvalidAllele { seq_id: String, allele: String },

    #[error("invalid strand '{0}': must be '+'/'1' or '-'/' -1'")]
    InvalidStrand(String),

    #[error("required field '{0}' cannot be empty or zero")]
    EmptyOrZeroField(String),

    #[error("malformed numeric field '{0}'")]
    MalformedNumber(String),

    #[error("record kind {0:?} is evidence or non-product record, not a product mutation")]
    NotProductMutation(GdKind),

    #[error("invalid DV1 EventId format: '{0}' (must match ^DV1_[0-9a-f]{{32}}$)")]
    InvalidEventId(String),

    #[error("hash collision detected for EventId '{event_id}' with different preimages")]
    EventIdCollision {
        event_id: String,
        first_preimage: Vec<u8>,
        second_preimage: Vec<u8>,
    },

    #[error("differential construction invariant failed: {0}")]
    InvariantViolation(&'static str),
}

/// Content-derived stable event identifier (DV1_<32 lowercase hex chars>).
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EventId(String);

impl EventId {
    pub fn parse(s: &str) -> Result<Self, DifferentialError> {
        if s.len() != 36 || !s.starts_with("DV1_") {
            return Err(DifferentialError::InvalidEventId(s.to_string()));
        }
        let hex_part = &s[4..];
        for b in hex_part.bytes() {
            if !b.is_ascii_hexdigit() || b.is_ascii_uppercase() {
                return Err(DifferentialError::InvalidEventId(s.to_string()));
            }
        }
        Ok(Self(s.to_string()))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for EventId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// Canonical reference contig identified by exact seq_id and 32-byte sequence SHA-256.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CanonicalContig {
    pub seq_id: String,
    pub sequence_sha256: [u8; 32],
}

impl CanonicalContig {
    pub fn from_ref_contig(ref_contig: &RefContig) -> Result<Self, DifferentialError> {
        if ref_contig.seq.is_empty() {
            return Err(DifferentialError::InvalidContigSequence {
                seq_id: ref_contig.name.clone(),
                msg: "sequence is empty".to_string(),
            });
        }
        let mut hasher = Sha256::new();
        for &b in &ref_contig.seq {
            if !b.is_ascii_alphabetic() && b != b'*' && b != b'-' {
                return Err(DifferentialError::InvalidContigSequence {
                    seq_id: ref_contig.name.clone(),
                    msg: format!("invalid sequence byte: {b}"),
                });
            }
            hasher.update([b.to_ascii_uppercase()]);
        }
        let digest: [u8; 32] = hasher.finalize().into();
        Ok(Self {
            seq_id: ref_contig.name.clone(),
            sequence_sha256: digest,
        })
    }
}

/// Canonical side of a sequence junction.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CanonicalJunctionSide {
    pub contig: CanonicalContig,
    pub position: u64,
    pub strand: i8, // exactly 1 or -1
}

/// Canonical physical event representation.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum CanonicalEvent {
    Snp {
        contig: CanonicalContig,
        position: u64,
        new_base: u8,
    },
    Sub {
        contig: CanonicalContig,
        position: u64,
        size: u64,
        new_seq: Vec<u8>,
    },
    Ins {
        contig: CanonicalContig,
        position: u64,
        new_seq: Vec<u8>,
    },
    Del {
        contig: CanonicalContig,
        position: u64,
        size: u64,
    },
    Mob {
        contig: CanonicalContig,
        position: u64,
        repeat_name: String,
        strand: i8,
        duplication_size: i64,
    },
    Amp {
        contig: CanonicalContig,
        position: u64,
        size: u64,
        new_copy_number: u64,
    },
    Con {
        contig: CanonicalContig,
        position: u64,
        size: u64,
        region: String,
    },
    Inv {
        contig: CanonicalContig,
        position: u64,
        size: u64,
    },
    Jc {
        side_1: CanonicalJunctionSide,
        side_2: CanonicalJunctionSide,
        overlap: i64,
    },
}

impl CanonicalEvent {
    /// Domain separation prefix for DV1 hash preimage.
    pub const DOMAIN_PREFIX: &'static [u8] = b"ProkaDiff\0DV1\0";

    /// Encode canonical event into its exact deterministic wire format.
    pub fn to_preimage_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(128);
        buf.extend_from_slice(Self::DOMAIN_PREFIX);

        match self {
            Self::Snp {
                contig,
                position,
                new_base,
            } => {
                buf.push(0x01);
                encode_contig(&mut buf, contig);
                encode_u64(&mut buf, *position);
                buf.push(*new_base);
            }
            Self::Sub {
                contig,
                position,
                size,
                new_seq,
            } => {
                buf.push(0x02);
                encode_contig(&mut buf, contig);
                encode_u64(&mut buf, *position);
                encode_u64(&mut buf, *size);
                encode_bytes(&mut buf, new_seq);
            }
            Self::Ins {
                contig,
                position,
                new_seq,
            } => {
                buf.push(0x03);
                encode_contig(&mut buf, contig);
                encode_u64(&mut buf, *position);
                encode_bytes(&mut buf, new_seq);
            }
            Self::Del {
                contig,
                position,
                size,
            } => {
                buf.push(0x04);
                encode_contig(&mut buf, contig);
                encode_u64(&mut buf, *position);
                encode_u64(&mut buf, *size);
            }
            Self::Mob {
                contig,
                position,
                repeat_name,
                strand,
                duplication_size,
            } => {
                buf.push(0x05);
                encode_contig(&mut buf, contig);
                encode_u64(&mut buf, *position);
                encode_bytes(&mut buf, repeat_name.as_bytes());
                encode_strand(&mut buf, *strand);
                encode_i64(&mut buf, *duplication_size);
            }
            Self::Amp {
                contig,
                position,
                size,
                new_copy_number,
            } => {
                buf.push(0x06);
                encode_contig(&mut buf, contig);
                encode_u64(&mut buf, *position);
                encode_u64(&mut buf, *size);
                encode_u64(&mut buf, *new_copy_number);
            }
            Self::Con {
                contig,
                position,
                size,
                region,
            } => {
                buf.push(0x07);
                encode_contig(&mut buf, contig);
                encode_u64(&mut buf, *position);
                encode_u64(&mut buf, *size);
                encode_bytes(&mut buf, region.as_bytes());
            }
            Self::Inv {
                contig,
                position,
                size,
            } => {
                buf.push(0x08);
                encode_contig(&mut buf, contig);
                encode_u64(&mut buf, *position);
                encode_u64(&mut buf, *size);
            }
            Self::Jc {
                side_1,
                side_2,
                overlap,
            } => {
                buf.push(0x09);
                encode_side(&mut buf, side_1);
                encode_side(&mut buf, side_2);
                encode_i64(&mut buf, *overlap);
            }
        }
        buf
    }

    /// Compute the content-derived DV1 EventId and return (id, preimage_bytes).
    pub fn compute_event_id(&self) -> (EventId, Vec<u8>) {
        let preimage = self.to_preimage_bytes();
        let digest = Sha256::digest(&preimage);
        let mut hex_str = String::with_capacity(36);
        hex_str.push_str("DV1_");
        const HEX: &[u8; 16] = b"0123456789abcdef";
        for &b in &digest[..16] {
            hex_str.push(char::from(HEX[usize::from(b >> 4)]));
            hex_str.push(char::from(HEX[usize::from(b & 0x0f)]));
        }
        (EventId(hex_str), preimage)
    }

    /// Construct a CanonicalEvent from a GenomeDiff entry and reference contig skeleton.
    pub fn from_gd_entry(entry: &GdEntry, refs: &[RefContig]) -> Result<Self, DifferentialError> {
        match entry.kind {
            GdKind::Snp => {
                let seq_id = entry
                    .fields
                    .first()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("seq_id".to_string()))?;
                let ref_contig = find_ref(refs, seq_id)?;
                let canonical_contig = CanonicalContig::from_ref_contig(ref_contig)?;
                let pos: u64 = parse_u64(entry.fields.get(1), "position")?;
                validate_position(seq_id, pos, ref_contig.seq.len())?;
                let allele = entry
                    .fields
                    .get(2)
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("allele".to_string()))?;
                let normalized = normalize_allele(seq_id, allele)?;
                let [b] = normalized.as_slice() else {
                    return Err(DifferentialError::InvalidAllele {
                        seq_id: seq_id.clone(),
                        allele: allele.clone(),
                    });
                };
                Ok(Self::Snp {
                    contig: canonical_contig,
                    position: pos,
                    new_base: *b,
                })
            }
            GdKind::Sub => {
                let seq_id = entry
                    .fields
                    .first()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("seq_id".to_string()))?;
                let ref_contig = find_ref(refs, seq_id)?;
                let canonical_contig = CanonicalContig::from_ref_contig(ref_contig)?;
                let pos: u64 = parse_u64(entry.fields.get(1), "position")?;
                let size: u64 = parse_u64(entry.fields.get(2), "size")?;
                validate_span(seq_id, pos, size, ref_contig.seq.len())?;
                let new_seq_str = entry
                    .fields
                    .get(3)
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("new_seq".to_string()))?;
                let new_seq = normalize_allele(seq_id, new_seq_str)?;
                Ok(Self::Sub {
                    contig: canonical_contig,
                    position: pos,
                    size,
                    new_seq,
                })
            }
            GdKind::Ins => {
                let seq_id = entry
                    .fields
                    .first()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("seq_id".to_string()))?;
                let ref_contig = find_ref(refs, seq_id)?;
                let canonical_contig = CanonicalContig::from_ref_contig(ref_contig)?;
                let pos: u64 = parse_u64(entry.fields.get(1), "position")?;
                validate_position(seq_id, pos, ref_contig.seq.len())?;
                let ins_seq_str = entry
                    .fields
                    .get(2)
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("new_seq".to_string()))?;
                let ins_raw = normalize_allele(seq_id, ins_seq_str)?;
                let (norm_pos, norm_seq) = right_align_ins(&ref_contig.seq, pos, &ins_raw);
                Ok(Self::Ins {
                    contig: canonical_contig,
                    position: norm_pos,
                    new_seq: norm_seq,
                })
            }
            GdKind::Del => {
                let seq_id = entry
                    .fields
                    .first()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("seq_id".to_string()))?;
                let ref_contig = find_ref(refs, seq_id)?;
                let canonical_contig = CanonicalContig::from_ref_contig(ref_contig)?;
                let pos: u64 = parse_u64(entry.fields.get(1), "position")?;
                let size: u64 = parse_u64(entry.fields.get(2), "size")?;
                validate_span(seq_id, pos, size, ref_contig.seq.len())?;
                let norm_pos = right_align_del(&ref_contig.seq, pos, size);
                Ok(Self::Del {
                    contig: canonical_contig,
                    position: norm_pos,
                    size,
                })
            }
            GdKind::Mob => {
                let seq_id = entry
                    .fields
                    .first()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("seq_id".to_string()))?;
                let ref_contig = find_ref(refs, seq_id)?;
                let canonical_contig = CanonicalContig::from_ref_contig(ref_contig)?;
                let pos: u64 = parse_u64(entry.fields.get(1), "position")?;
                validate_position(seq_id, pos, ref_contig.seq.len())?;
                let repeat_name = entry.fields.get(2).cloned().ok_or_else(|| {
                    DifferentialError::EmptyOrZeroField("repeat_name".to_string())
                })?;
                if repeat_name.is_empty() {
                    return Err(DifferentialError::EmptyOrZeroField(
                        "repeat_name".to_string(),
                    ));
                }
                let strand_token = entry
                    .fields
                    .get(3)
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("strand".to_string()))?;
                let strand = parse_strand(strand_token)?;
                let duplication_size: i64 = entry
                    .fields
                    .get(4)
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| {
                        DifferentialError::MalformedNumber("duplication_size".to_string())
                    })?;
                Ok(Self::Mob {
                    contig: canonical_contig,
                    position: pos,
                    repeat_name,
                    strand,
                    duplication_size,
                })
            }
            GdKind::Amp => {
                let seq_id = entry
                    .fields
                    .first()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("seq_id".to_string()))?;
                let ref_contig = find_ref(refs, seq_id)?;
                let canonical_contig = CanonicalContig::from_ref_contig(ref_contig)?;
                let pos: u64 = parse_u64(entry.fields.get(1), "position")?;
                let size: u64 = parse_u64(entry.fields.get(2), "size")?;
                let new_copy_number: u64 = parse_u64(entry.fields.get(3), "new_copy_number")?;
                validate_span(seq_id, pos, size, ref_contig.seq.len())?;
                if new_copy_number == 0 {
                    return Err(DifferentialError::EmptyOrZeroField(
                        "new_copy_number".to_string(),
                    ));
                }
                Ok(Self::Amp {
                    contig: canonical_contig,
                    position: pos,
                    size,
                    new_copy_number,
                })
            }
            GdKind::Con => {
                let seq_id = entry
                    .fields
                    .first()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("seq_id".to_string()))?;
                let ref_contig = find_ref(refs, seq_id)?;
                let canonical_contig = CanonicalContig::from_ref_contig(ref_contig)?;
                let pos: u64 = parse_u64(entry.fields.get(1), "position")?;
                let size: u64 = parse_u64(entry.fields.get(2), "size")?;
                let region = entry
                    .fields
                    .get(3)
                    .cloned()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("region".to_string()))?;
                validate_span(seq_id, pos, size, ref_contig.seq.len())?;
                if region.is_empty() {
                    return Err(DifferentialError::EmptyOrZeroField("region".to_string()));
                }
                Ok(Self::Con {
                    contig: canonical_contig,
                    position: pos,
                    size,
                    region,
                })
            }
            GdKind::Inv => {
                let seq_id = entry
                    .fields
                    .first()
                    .ok_or_else(|| DifferentialError::EmptyOrZeroField("seq_id".to_string()))?;
                let ref_contig = find_ref(refs, seq_id)?;
                let canonical_contig = CanonicalContig::from_ref_contig(ref_contig)?;
                let pos: u64 = parse_u64(entry.fields.get(1), "position")?;
                let size: u64 = parse_u64(entry.fields.get(2), "size")?;
                validate_span(seq_id, pos, size, ref_contig.seq.len())?;
                Ok(Self::Inv {
                    contig: canonical_contig,
                    position: pos,
                    size,
                })
            }
            GdKind::Jc => {
                let s1_id = entry.fields.first().ok_or_else(|| {
                    DifferentialError::EmptyOrZeroField("side1 seq_id".to_string())
                })?;
                let ref_1 = find_ref(refs, s1_id)?;
                let contig_1 = CanonicalContig::from_ref_contig(ref_1)?;
                let pos_1: u64 = parse_u64(entry.fields.get(1), "side1 position")?;
                let strand_1 = parse_strand(entry.fields.get(2).map_or("", String::as_str))?;
                validate_position(s1_id, pos_1, ref_1.seq.len())?;

                let s2_id = entry.fields.get(3).ok_or_else(|| {
                    DifferentialError::EmptyOrZeroField("side2 seq_id".to_string())
                })?;
                let ref_2 = find_ref(refs, s2_id)?;
                let contig_2 = CanonicalContig::from_ref_contig(ref_2)?;
                let pos_2: u64 = parse_u64(entry.fields.get(4), "side2 position")?;
                let strand_2 = parse_strand(entry.fields.get(5).map_or("", String::as_str))?;
                validate_position(s2_id, pos_2, ref_2.seq.len())?;

                let overlap: i64 = entry
                    .fields
                    .get(6)
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| DifferentialError::MalformedNumber("overlap".to_string()))?;

                let side_a = CanonicalJunctionSide {
                    contig: contig_1,
                    position: pos_1,
                    strand: strand_1,
                };
                let side_b = CanonicalJunctionSide {
                    contig: contig_2,
                    position: pos_2,
                    strand: strand_2,
                };

                let (side_1, side_2) = if side_a <= side_b {
                    (side_a, side_b)
                } else {
                    (side_b, side_a)
                };

                Ok(Self::Jc {
                    side_1,
                    side_2,
                    overlap,
                })
            }
            other => Err(DifferentialError::NotProductMutation(other)),
        }
    }
}

/// Collision-guard for registering computed event IDs during construction.
#[derive(Default, Debug, Clone)]
pub struct EventIdCollisionGuard {
    seen: BTreeMap<EventId, Vec<u8>>,
}

impl EventIdCollisionGuard {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn check_or_insert(
        &mut self,
        id: EventId,
        preimage: Vec<u8>,
    ) -> Result<EventId, DifferentialError> {
        if let Some(existing) = self.seen.get(&id) {
            if existing != &preimage {
                return Err(DifferentialError::EventIdCollision {
                    event_id: id.0,
                    first_preimage: existing.clone(),
                    second_preimage: preimage,
                });
            }
        } else {
            self.seen.insert(id.clone(), preimage);
        }
        Ok(id)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum EvidenceKind {
    Ra,
    Mc,
    Jc,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct GdEvidenceRef {
    pub gd_id: u32,
    pub kind: EvidenceKind,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EvidenceReferences {
    pub gd_parent_ids: Option<Vec<u32>>,
    pub ra: Option<Vec<GdEvidenceRef>>,
    pub mc: Option<Vec<GdEvidenceRef>>,
    pub jc: Option<Vec<GdEvidenceRef>>,
    pub pd_attributes: BTreeMap<String, Vec<String>>,
    pub unresolved_parent_ids: Vec<u32>,
}

impl EvidenceReferences {
    /// Resolve evidence references and `pd_*` attributes for a group of GD entries.
    ///
    /// - `group_entries`: The GD entries merged into this canonical event.
    /// - `all_records`: Optional slice of all GD records in the document (for looking up parent records).
    pub fn resolve_for_group(group_entries: &[&GdEntry], all_records: Option<&[GdEntry]>) -> Self {
        let mut all_parents: Vec<u32> = group_entries
            .iter()
            .flat_map(|e| e.parent_ids.iter().copied())
            .collect();
        all_parents.sort_unstable();
        all_parents.dedup();

        let mut ra_refs: Vec<GdEvidenceRef> = Vec::new();
        let mut mc_refs: Vec<GdEvidenceRef> = Vec::new();
        let mut jc_refs: Vec<GdEvidenceRef> = Vec::new();
        let mut unresolved: Vec<u32> = Vec::new();

        let mut pd_attrs: BTreeMap<String, Vec<String>> = BTreeMap::new();

        // 1. Collect pd_* attributes from the mutation group entries themselves.
        for e in group_entries {
            for (k, v) in &e.attrs {
                if k.starts_with("pd_") {
                    pd_attrs.entry(k.clone()).or_default().push(v.clone());
                }
            }
        }

        // 2. A differential JC may reference its own actual JC record.
        for e in group_entries {
            if e.kind == GdKind::Jc {
                jc_refs.push(GdEvidenceRef {
                    gd_id: e.id,
                    kind: EvidenceKind::Jc,
                });
            }
        }

        // 3. Resolve parent_ids against all_records.
        for &pid in &all_parents {
            if let Some(parent) = all_records.and_then(|recs| recs.iter().find(|r| r.id == pid)) {
                match parent.kind {
                    GdKind::Ra => {
                        ra_refs.push(GdEvidenceRef {
                            gd_id: pid,
                            kind: EvidenceKind::Ra,
                        });
                    }
                    GdKind::Mc => {
                        mc_refs.push(GdEvidenceRef {
                            gd_id: pid,
                            kind: EvidenceKind::Mc,
                        });
                    }
                    GdKind::Jc => {
                        jc_refs.push(GdEvidenceRef {
                            gd_id: pid,
                            kind: EvidenceKind::Jc,
                        });
                    }
                    _ => {
                        unresolved.push(pid);
                    }
                }
                // Also collect pd_* attributes from resolved evidence records.
                for (k, v) in &parent.attrs {
                    if k.starts_with("pd_") {
                        pd_attrs.entry(k.clone()).or_default().push(v.clone());
                    }
                }
            } else {
                unresolved.push(pid);
            }
        }

        ra_refs.sort_unstable();
        ra_refs.dedup();
        mc_refs.sort_unstable();
        mc_refs.dedup();
        jc_refs.sort_unstable();
        jc_refs.dedup();
        unresolved.sort_unstable();
        unresolved.dedup();

        for vals in pd_attrs.values_mut() {
            vals.sort_unstable();
            vals.dedup();
        }

        Self {
            gd_parent_ids: if all_parents.is_empty() {
                None
            } else {
                Some(all_parents)
            },
            ra: if ra_refs.is_empty() {
                None
            } else {
                Some(ra_refs)
            },
            mc: if mc_refs.is_empty() {
                None
            } else {
                Some(mc_refs)
            },
            jc: if jc_refs.is_empty() {
                None
            } else {
                Some(jc_refs)
            },
            pd_attributes: pd_attrs,
            unresolved_parent_ids: unresolved,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DifferentialEvent {
    pub event_id: EventId,
    pub canonical: CanonicalEvent,
    pub representative: GdEntry,
    pub merged_source_ids: Vec<u32>,
    pub evidence: EvidenceReferences,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DifferentialResultSet {
    pub events: Vec<DifferentialEvent>,
    pub starter_mutation_count: usize,
}

#[derive(Clone, Debug)]
struct CanonicalGroup {
    canonical: CanonicalEvent,
    canonical_bytes: Vec<u8>,
    representative: GdEntry,
    merged_source_ids: Vec<u32>,
    evidence: EvidenceReferences,
}

type RepresentativeSortKey = (
    u8,
    Vec<u8>,
    Vec<String>,
    u32,
    Vec<u32>,
    Vec<(String, String)>,
);

fn representative_sort_key(
    entry: &GdEntry,
    kind_tag: u8,
    canonical_bytes: &[u8],
) -> RepresentativeSortKey {
    let mut parent_ids = entry.parent_ids.clone();
    parent_ids.sort_unstable();
    parent_ids.dedup();
    let attrs: Vec<(String, String)> = entry
        .attrs
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    (
        kind_tag,
        canonical_bytes.to_vec(),
        entry.fields.clone(),
        entry.id,
        parent_ids,
        attrs,
    )
}

fn coalesce_canonical_entries(
    entries: &[&GdEntry],
    refs: &[RefContig],
    all_records: Option<&[GdEntry]>,
) -> Result<Vec<CanonicalGroup>, DifferentialError> {
    let mut map: BTreeMap<Vec<u8>, (CanonicalEvent, Vec<GdEntry>)> = BTreeMap::new();
    for entry in entries {
        let canonical = CanonicalEvent::from_gd_entry(entry, refs)?;
        let bytes = canonical.to_preimage_bytes();
        map.entry(bytes)
            .or_insert_with(|| (canonical, Vec::new()))
            .1
            .push((*entry).clone());
    }

    let mut groups = Vec::with_capacity(map.len());
    for (canonical_bytes, (canonical, group_entries)) in map {
        let kind_tag = canonical_bytes
            .get(CanonicalEvent::DOMAIN_PREFIX.len())
            .copied()
            .ok_or(DifferentialError::InvariantViolation(
                "canonical preimage is missing its kind tag",
            ))?;
        let rep = group_entries
            .iter()
            .min_by_key(|e| representative_sort_key(e, kind_tag, &canonical_bytes))
            .ok_or(DifferentialError::InvariantViolation(
                "canonical group has no source entry",
            ))?
            .clone();

        let mut merged_source_ids: Vec<u32> = group_entries.iter().map(|e| e.id).collect();
        merged_source_ids.sort_unstable();
        merged_source_ids.dedup();

        let group_refs: Vec<&GdEntry> = group_entries.iter().collect();
        let evidence = EvidenceReferences::resolve_for_group(&group_refs, all_records);

        groups.push(CanonicalGroup {
            canonical,
            canonical_bytes,
            representative: rep,
            merged_source_ids,
            evidence,
        });
    }

    groups.sort_by(|a, b| a.canonical_bytes.cmp(&b.canonical_bytes));
    Ok(groups)
}

/// Compute the canonical differential event set: M_edited \ M_starter.
pub fn build_differential_events(
    starter_gd: &GenomeDiff,
    edited_gd: &GenomeDiff,
    refs: &[RefContig],
) -> Result<DifferentialResultSet, DifferentialError> {
    let starter_absorbed_jc_ids = mob_constituent_jc_ids(starter_gd);
    let edited_absorbed_jc_ids = mob_constituent_jc_ids(edited_gd);

    let starter_product: Vec<&GdEntry> = starter_gd
        .entries
        .iter()
        .filter(|e| is_product_mutation(e.kind) && !starter_absorbed_jc_ids.contains(&e.id))
        .collect();

    let edited_product: Vec<&GdEntry> = edited_gd
        .entries
        .iter()
        .filter(|e| is_product_mutation(e.kind) && !edited_absorbed_jc_ids.contains(&e.id))
        .collect();

    let starter_groups =
        coalesce_canonical_entries(&starter_product, refs, Some(&starter_gd.entries))?;
    let edited_groups =
        coalesce_canonical_entries(&edited_product, refs, Some(&edited_gd.entries))?;

    let starter_mutation_count = starter_groups.len();

    let mut used_starter = vec![false; starter_groups.len()];
    let mut subtracted_edited = vec![false; edited_groups.len()];

    // Phase 1: One-to-one exact canonical match subtraction
    for (e_idx, edited) in edited_groups.iter().enumerate() {
        for (s_idx, starter) in starter_groups.iter().enumerate() {
            if !used_starter[s_idx] && starter.canonical_bytes == edited.canonical_bytes {
                used_starter[s_idx] = true;
                subtracted_edited[e_idx] = true;
                break;
            }
        }
    }

    let mut candidates = Vec::new();
    for (edited_index, edited) in edited_groups.iter().enumerate() {
        if subtracted_edited[edited_index] {
            continue;
        }
        for (starter_index, starter) in starter_groups.iter().enumerate() {
            if used_starter[starter_index] {
                continue;
            }
            if let Some(cost) = tolerant_match_cost(&edited.canonical, &starter.canonical) {
                candidates.push(TolerantCandidate {
                    edited_index,
                    starter_index,
                    cost,
                });
            }
        }
    }
    candidates.sort_unstable();
    for (edited_index, starter_index) in
        minimum_cost_maximum_matching(edited_groups.len(), starter_groups.len(), &candidates)
    {
        if !subtracted_edited[edited_index] && !used_starter[starter_index] {
            subtracted_edited[edited_index] = true;
            used_starter[starter_index] = true;
        }
    }

    let mut guard = EventIdCollisionGuard::new();
    let mut events = Vec::new();
    for (e_idx, edited) in edited_groups.into_iter().enumerate() {
        if !subtracted_edited[e_idx] {
            let (event_id, preimage) = edited.canonical.compute_event_id();
            guard.check_or_insert(event_id.clone(), preimage)?;
            events.push(DifferentialEvent {
                event_id,
                canonical: edited.canonical,
                representative: edited.representative,
                merged_source_ids: edited.merged_source_ids,
                evidence: edited.evidence,
            });
        }
    }

    events.sort_by(|a, b| a.canonical.cmp(&b.canonical));

    Ok(DifferentialResultSet {
        events,
        starter_mutation_count,
    })
}

fn mob_constituent_jc_ids(gd: &GenomeDiff) -> HashSet<u32> {
    gd.entries
        .iter()
        .filter(|entry| entry.kind == GdKind::Mob)
        .flat_map(|entry| entry.parent_ids.iter().copied())
        .filter(|parent_id| {
            gd.entries
                .iter()
                .any(|candidate| candidate.id == *parent_id && candidate.kind == GdKind::Jc)
        })
        .collect()
}

fn tolerant_match_cost(
    edited: &CanonicalEvent,
    starter: &CanonicalEvent,
) -> Option<CoordinateCost> {
    let deltas = match (edited, starter) {
        (
            CanonicalEvent::Jc {
                side_1: edited_1,
                side_2: edited_2,
                overlap: edited_overlap,
            },
            CanonicalEvent::Jc {
                side_1: starter_1,
                side_2: starter_2,
                overlap: starter_overlap,
            },
        ) if edited_overlap == starter_overlap
            && edited_1.contig == starter_1.contig
            && edited_1.strand == starter_1.strand
            && edited_2.contig == starter_2.contig
            && edited_2.strand == starter_2.strand =>
        {
            [
                edited_1.position.abs_diff(starter_1.position),
                edited_2.position.abs_diff(starter_2.position),
            ]
        }
        (
            CanonicalEvent::Mob {
                contig: edited_contig,
                position: edited_position,
                repeat_name: edited_repeat,
                strand: edited_strand,
                duplication_size: edited_duplication,
            },
            CanonicalEvent::Mob {
                contig: starter_contig,
                position: starter_position,
                repeat_name: starter_repeat,
                strand: starter_strand,
                duplication_size: starter_duplication,
            },
        ) if edited_contig == starter_contig
            && edited_repeat == starter_repeat
            && edited_strand == starter_strand
            && edited_duplication == starter_duplication =>
        {
            [
                edited_position.abs_diff(*starter_position),
                edited_position.abs_diff(*starter_position),
            ]
        }
        (
            CanonicalEvent::Del {
                contig: edited_contig,
                position: edited_position,
                size: edited_size,
            },
            CanonicalEvent::Del {
                contig: starter_contig,
                position: starter_position,
                size: starter_size,
            },
        ) if edited_contig == starter_contig && *edited_size > 2 && *starter_size > 2 => {
            let edited_end = interval_end(*edited_position, *edited_size)?;
            let starter_end = interval_end(*starter_position, *starter_size)?;
            [
                edited_position.abs_diff(*starter_position),
                edited_end.abs_diff(starter_end),
            ]
        }
        _ => return None,
    };
    if deltas.iter().any(|delta| *delta > STRUCTURAL_TOLERANCE_BP) {
        return None;
    }
    let max_delta = *deltas.iter().max()?;
    let sum_delta = deltas[0].checked_add(deltas[1])?;
    Some(CoordinateCost {
        max_delta: i64::try_from(max_delta).ok()?,
        sum_delta: i64::try_from(sum_delta).ok()?,
    })
}

fn interval_end(position: u64, size: u64) -> Option<u64> {
    position.checked_add(size.checked_sub(1)?)
}

fn encode_contig(buf: &mut Vec<u8>, contig: &CanonicalContig) {
    let bytes = contig.seq_id.as_bytes();
    buf.extend_from_slice(&(bytes.len() as u32).to_be_bytes());
    buf.extend_from_slice(bytes);
    buf.extend_from_slice(&contig.sequence_sha256);
}

fn encode_u64(buf: &mut Vec<u8>, val: u64) {
    buf.extend_from_slice(&val.to_be_bytes());
}

fn encode_i64(buf: &mut Vec<u8>, val: i64) {
    buf.extend_from_slice(&val.to_be_bytes());
}

fn encode_bytes(buf: &mut Vec<u8>, bytes: &[u8]) {
    buf.extend_from_slice(&(bytes.len() as u32).to_be_bytes());
    buf.extend_from_slice(bytes);
}

fn encode_strand(buf: &mut Vec<u8>, strand: i8) {
    if strand == 1 {
        buf.push(0x01);
    } else {
        buf.push(0x02);
    }
}

fn encode_side(buf: &mut Vec<u8>, side: &CanonicalJunctionSide) {
    encode_contig(buf, &side.contig);
    encode_u64(buf, side.position);
    encode_strand(buf, side.strand);
}

fn find_ref<'a>(refs: &'a [RefContig], seq_id: &str) -> Result<&'a RefContig, DifferentialError> {
    refs.iter()
        .find(|r| r.name == seq_id)
        .ok_or_else(|| DifferentialError::ContigNotFound(seq_id.to_string()))
}

fn validate_position(
    seq_id: &str,
    position: u64,
    contig_len: usize,
) -> Result<(), DifferentialError> {
    let contig_len_u64 =
        u64::try_from(contig_len).map_err(|_| DifferentialError::InvalidContigSequence {
            seq_id: seq_id.to_string(),
            msg: "sequence length cannot be represented as a genomic coordinate".to_string(),
        })?;
    if position == 0 || position > contig_len_u64 {
        return Err(DifferentialError::InvalidPosition {
            seq_id: seq_id.to_string(),
            position,
            contig_len,
        });
    }
    Ok(())
}

fn validate_span(
    seq_id: &str,
    position: u64,
    size: u64,
    contig_len: usize,
) -> Result<(), DifferentialError> {
    validate_position(seq_id, position, contig_len)?;
    let contig_len_u64 =
        u64::try_from(contig_len).map_err(|_| DifferentialError::InvalidContigSequence {
            seq_id: seq_id.to_string(),
            msg: "sequence length cannot be represented as a genomic coordinate".to_string(),
        })?;
    let end_exclusive = position
        .checked_sub(1)
        .and_then(|start| start.checked_add(size));
    if size == 0 || end_exclusive.is_none_or(|end| end > contig_len_u64) {
        return Err(DifferentialError::InvalidSize {
            seq_id: seq_id.to_string(),
            position,
            size,
            contig_len,
        });
    }
    Ok(())
}

fn normalize_allele(seq_id: &str, allele: &str) -> Result<Vec<u8>, DifferentialError> {
    if allele.is_empty() {
        return Err(DifferentialError::InvalidAllele {
            seq_id: seq_id.to_string(),
            allele: allele.to_string(),
        });
    }
    let normalized: Vec<u8> = allele
        .bytes()
        .map(|base| base.to_ascii_uppercase())
        .collect();
    if normalized
        .iter()
        .any(|base| !DNA_IUPAC_ALPHABET.contains(base))
    {
        return Err(DifferentialError::InvalidAllele {
            seq_id: seq_id.to_string(),
            allele: allele.to_string(),
        });
    }
    Ok(normalized)
}

fn parse_u64(s: Option<&String>, field: &str) -> Result<u64, DifferentialError> {
    s.and_then(|x| x.parse().ok())
        .ok_or_else(|| DifferentialError::MalformedNumber(field.to_string()))
}

fn parse_strand(s: &str) -> Result<i8, DifferentialError> {
    match s {
        "+" | "1" => Ok(1),
        "-" | "-1" => Ok(-1),
        _ => Err(DifferentialError::InvalidStrand(s.to_string())),
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    fn test_ref_skeleton() -> Vec<RefContig> {
        vec![
            RefContig {
                name: "chr1".into(),
                seq: b"GCTGGTAAACCATAACTGTCGCAGCGTTAACCGGTT".to_vec(),
            },
            RefContig {
                name: "chr2".into(),
                seq: b"ATCGATCGATCGATCGATCG".to_vec(),
            },
        ]
    }

    #[test]
    fn dv1_canonical_golden_hashes() {
        let refs = test_ref_skeleton();

        // 1. SNP
        let snp_entry = GdEntry::snp(1, "chr1", 5, "C");
        let snp_event = CanonicalEvent::from_gd_entry(&snp_entry, &refs).unwrap();
        let (snp_id, snp_preimage) = snp_event.compute_event_id();
        assert_eq!(snp_id.as_str(), "DV1_1ceede1cac40a9b5c09cb37ca592b1a2");
        assert_eq!(snp_id.as_str().len(), 36);
        assert_eq!(&snp_preimage[..14], b"ProkaDiff\0DV1\0");
        assert_eq!(snp_preimage[14], 0x01); // kind SNP

        // 2. INS with 3' right-alignment
        let ins_entry = GdEntry::ins(2, "chr1", 10, "A");
        let ins_event = CanonicalEvent::from_gd_entry(&ins_entry, &refs).unwrap();
        let (ins_id, _) = ins_event.compute_event_id();
        assert_eq!(ins_id.as_str(), "DV1_537bbf2a14670ce5f104d85dbd97708d");

        // 3. DEL with 3' right-alignment
        let del_entry = GdEntry::del(3, "chr1", 5, 2);
        let del_event = CanonicalEvent::from_gd_entry(&del_entry, &refs).unwrap();
        let (del_id, _) = del_event.compute_event_id();
        assert_eq!(del_id.as_str(), "DV1_9fb25881125d4c8207e7a6a34725640c");

        // 4. MOB
        let mob_entry = GdEntry::mob(4, "chr1", 8, "IS1", "+", 4);
        let mob_event = CanonicalEvent::from_gd_entry(&mob_entry, &refs).unwrap();
        let (mob_id, _) = mob_event.compute_event_id();
        assert_eq!(mob_id.as_str(), "DV1_ce91da4e6b3ae24faf1fb507f300a1c7");

        // 5. Reversed JC sides converge to identical event_id
        let jc1 = GdEntry::jc(5, "chr1", 10, "+", "chr2", 5, "-", 0);
        let jc2 = GdEntry::jc(6, "chr2", 5, "-", "chr1", 10, "+", 0);
        let event_jc1 = CanonicalEvent::from_gd_entry(&jc1, &refs).unwrap();
        let event_jc2 = CanonicalEvent::from_gd_entry(&jc2, &refs).unwrap();
        let (id1, pre1) = event_jc1.compute_event_id();
        let (id2, pre2) = event_jc2.compute_event_id();
        assert_eq!(id1.as_str(), "DV1_af8969c38915f6797f4204c6a93b08b4");
        assert_eq!(
            id1, id2,
            "reversed JC sides must yield identical DV1 EventId"
        );
        assert_eq!(
            pre1, pre2,
            "reversed JC sides must yield identical wire bytes"
        );
    }

    #[test]
    fn dv1_excludes_nonidentity_fields() {
        let refs = test_ref_skeleton();
        let mut entry1 = GdEntry::snp(100, "chr1", 5, "C");
        entry1.parent_ids = vec![1, 2, 3];
        entry1.attrs.insert("frequency".into(), "0.99".into());
        entry1.attrs.insert("pd_ra_depth".into(), "42".into());

        let mut entry2 = GdEntry::snp(999, "chr1", 5, "C");
        entry2.parent_ids = vec![99];
        entry2.attrs.insert("note".into(), "unrelated".into());

        let event1 = CanonicalEvent::from_gd_entry(&entry1, &refs).unwrap();
        let event2 = CanonicalEvent::from_gd_entry(&entry2, &refs).unwrap();

        let (id1, pre1) = event1.compute_event_id();
        let (id2, pre2) = event2.compute_event_id();
        assert_eq!(id1, id2);
        assert_eq!(pre1, pre2);
    }

    #[test]
    fn dv1_rejects_invalid() {
        let refs = test_ref_skeleton();

        // Position 0
        let zero_pos = GdEntry::snp(1, "chr1", 0, "A");
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&zero_pos, &refs),
            Err(DifferentialError::InvalidPosition { .. })
        ));

        // Unknown contig
        let unk_contig = GdEntry::snp(1, "chrX", 5, "A");
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&unk_contig, &refs),
            Err(DifferentialError::ContigNotFound(_))
        ));

        // Out of range position
        let oob_pos = GdEntry::snp(1, "chr1", 9999, "A");
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&oob_pos, &refs),
            Err(DifferentialError::InvalidPosition { .. })
        ));

        // Size 0 DEL
        let zero_del = GdEntry::del(1, "chr1", 5, 0);
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&zero_del, &refs),
            Err(DifferentialError::InvalidSize { .. })
        ));

        // Invalid strand
        let bad_strand = GdEntry::mob(1, "chr1", 5, "IS1", "invalid", 0);
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&bad_strand, &refs),
            Err(DifferentialError::InvalidStrand(_))
        ));

        // Non-product record (RA)
        let ra_entry = GdEntry::ra(1, "chr1", 5, 0, "A", "C");
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&ra_entry, &refs),
            Err(DifferentialError::NotProductMutation(GdKind::Ra))
        ));
    }

    #[test]
    fn canonical_rejects_out_of_range() {
        let refs = test_ref_skeleton();
        let oob_pos = GdEntry::snp(1, "chr1", 9999, "A");
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&oob_pos, &refs),
            Err(DifferentialError::InvalidPosition { .. })
        ));
        let oob_size = GdEntry::del(1, "chr1", 30, 20);
        assert!(matches!(
            CanonicalEvent::from_gd_entry(&oob_size, &refs),
            Err(DifferentialError::InvalidSize { .. })
        ));
    }

    #[test]
    fn dv1_detects_collision() {
        let mut guard = EventIdCollisionGuard::new();
        let id = EventId("DV1_0123456789abcdef0123456789abcdef".into());
        let pre1 = b"preimage_one".to_vec();
        let pre2 = b"preimage_two".to_vec();

        assert!(guard.check_or_insert(id.clone(), pre1).is_ok());
        // Same ID with identical preimage is fine (idempotent)
        assert!(guard
            .check_or_insert(id.clone(), b"preimage_one".to_vec())
            .is_ok());
        // Same ID with different preimage triggers collision error
        let err = guard.check_or_insert(id, pre2);
        assert!(matches!(
            err,
            Err(DifferentialError::EventIdCollision { .. })
        ));
    }

    #[test]
    fn subtract() {
        let refs = test_ref_skeleton();

        let mut starter = GenomeDiff::new();
        // 1. Inherited SNP
        starter.entries.push(GdEntry::snp(1, "chr1", 5, "C"));
        // 2. Inherited JC
        starter
            .entries
            .push(GdEntry::jc(2, "chr1", 10, "+", "chr2", 5, "-", 0));
        // 3. Inherited MOB
        starter
            .entries
            .push(GdEntry::mob(3, "chr1", 8, "IS1", "+", 4));
        // 4. Inherited structural DEL
        starter.entries.push(GdEntry::del(4, "chr1", 15, 6));

        let mut edited = GenomeDiff::new();
        // 1. Exact SNP matching
        edited.entries.push(GdEntry::snp(10, "chr1", 5, "C"));
        // 2. Jittered JC (both sides jitter by <= 2 bp)
        edited
            .entries
            .push(GdEntry::jc(20, "chr1", 12, "+", "chr2", 7, "-", 0));
        // 3. Jittered MOB (position jitters by 2 bp)
        edited
            .entries
            .push(GdEntry::mob(30, "chr1", 10, "IS1", "+", 4));
        // 4. Jittered structural DEL (start and end jitter by 1 bp)
        edited.entries.push(GdEntry::del(40, "chr1", 16, 6));
        // 5. Genuine edited-only SNP
        edited.entries.push(GdEntry::snp(50, "chr1", 20, "A"));

        let res = build_differential_events(&starter, &edited, &refs).unwrap();
        assert_eq!(res.starter_mutation_count, 4);
        assert_eq!(res.events.len(), 1);
        assert_eq!(res.events[0].representative.id, 50);
        assert!(matches!(
            res.events[0].canonical,
            CanonicalEvent::Snp {
                position: 20,
                new_base: b'A',
                ..
            }
        ));
    }

    #[test]
    fn subtract_rejects_tolerance() {
        let refs = test_ref_skeleton();

        let mut starter = GenomeDiff::new();
        starter
            .entries
            .push(GdEntry::jc(1, "chr1", 10, "+", "chr2", 5, "-", 0));
        starter
            .entries
            .push(GdEntry::mob(2, "chr1", 8, "IS1", "+", 4));
        starter.entries.push(GdEntry::del(3, "chr1", 15, 6));
        starter.entries.push(GdEntry::del(4, "chr1", 1, 2)); // short DEL at 1
        starter.entries.push(GdEntry::snp(5, "chr1", 30, "A"));
        starter.entries.push(GdEntry::ins(6, "chr1", 1, "TT"));
        starter.entries.push(GdEntry::amp(7, "chr1", 1, 10, 2));
        starter.entries.push(GdEntry::inv(8, "chr1", 1, 10));

        let mut edited = GenomeDiff::new();
        // 1. JC beyond 5 bp (delta 6 on side 1)
        edited
            .entries
            .push(GdEntry::jc(10, "chr1", 16, "+", "chr2", 5, "-", 0));
        // 2. MOB beyond 5 bp (delta 6)
        edited
            .entries
            .push(GdEntry::mob(20, "chr1", 14, "IS1", "+", 4));
        // 3. Structural DEL beyond 5 bp (delta 6)
        edited.entries.push(GdEntry::del(30, "chr1", 21, 6));
        // 4. Short DEL within 5 bp (delta 1, but short DEL is exact-only!)
        edited.entries.push(GdEntry::del(40, "chr1", 2, 2));
        // 5. SNP with 1 bp delta -> exact only!
        edited.entries.push(GdEntry::snp(50, "chr1", 31, "A"));
        // 6. INS with 1 bp delta -> exact only!
        edited.entries.push(GdEntry::ins(60, "chr1", 2, "TT"));
        // 7. AMP with 1 bp delta -> exact only!
        edited.entries.push(GdEntry::amp(70, "chr1", 2, 10, 2));
        // 8. INV with 1 bp delta -> exact only!
        edited.entries.push(GdEntry::inv(80, "chr1", 2, 10));

        let res = build_differential_events(&starter, &edited, &refs).unwrap();
        // All 8 edited mutations must survive because tolerance is rejected for all of them!
        assert_eq!(res.events.len(), 8);
    }

    #[test]
    fn subtract_one_to_one() {
        let refs = test_ref_skeleton();

        let mut starter = GenomeDiff::new();
        // Single starter structural DEL
        starter.entries.push(GdEntry::del(1, "chr1", 15, 6));

        let mut edited = GenomeDiff::new();
        // Two candidate deletions close to starter:
        // pos 16 right-aligns to 17 (delta 2 from 15)
        // pos 18 right-aligns to 19 (delta 4 from 15)
        edited.entries.push(GdEntry::del(10, "chr1", 16, 6));
        edited.entries.push(GdEntry::del(20, "chr1", 18, 6));

        let res = build_differential_events(&starter, &edited, &refs).unwrap();
        assert_eq!(res.starter_mutation_count, 1);
        // Exactly one matches starter DEL, exactly one survives!
        assert_eq!(res.events.len(), 1);
        // DEL at 16 (delta 2) is closer and consumed, DEL at 18 (id 20) survives!
        assert_eq!(res.events[0].representative.id, 20);
    }

    #[test]
    fn subtract_nontransitive_chain() {
        let refs = test_ref_skeleton();

        let mut starter = GenomeDiff::new();
        starter.entries.push(GdEntry::del(1, "chr1", 15, 6));

        // Test with edited entries in order A then B
        let mut edited_ab = GenomeDiff::new();
        edited_ab.entries.push(GdEntry::del(10, "chr1", 12, 6)); // delta 3
        edited_ab.entries.push(GdEntry::del(20, "chr1", 18, 6)); // delta 3

        let res_ab = build_differential_events(&starter, &edited_ab, &refs).unwrap();

        // Test with edited entries in reversed order B then A
        let mut edited_ba = GenomeDiff::new();
        edited_ba.entries.push(GdEntry::del(20, "chr1", 18, 6)); // delta 3
        edited_ba.entries.push(GdEntry::del(10, "chr1", 12, 6)); // delta 3

        let res_ba = build_differential_events(&starter, &edited_ba, &refs).unwrap();

        assert_eq!(res_ab.events.len(), 1);
        assert_eq!(res_ba.events.len(), 1);
        // Regardless of input entry order, the resulting surviving event ID and content are identical!
        assert_eq!(res_ab.events[0].event_id, res_ba.events[0].event_id);
    }

    #[test]
    fn evidence() {
        use crate::audit::EvidenceSummary;

        let mut ra_entry = GdEntry::ra(10, "chr1", 5, 0, "A", "C");
        ra_entry.attrs.insert("pd_ra_depth".into(), "42.5".into());
        ra_entry
            .attrs
            .insert("pd_ra_support_reads".into(), "35".into());
        ra_entry
            .attrs
            .insert("pd_ra_frequency".into(), "0.95".into());

        let mut mc_entry = GdEntry::mc(20, "chr1", 10, 20, 10, 20);
        mc_entry.attrs.insert("pd_mc_cov".into(), "12".into());

        let mut jc_entry = GdEntry::jc(30, "chr1", 10, "+", "chr2", 5, "-", 0);
        jc_entry
            .attrs
            .insert("pd_support_reads".into(), "15".into());

        let mut mob_entry = GdEntry::mob(40, "chr1", 10, "IS1", "+", 4);
        mob_entry.parent_ids = vec![30]; // explicit parent JC

        let mut snp_entry = GdEntry::snp(50, "chr1", 5, "C");
        snp_entry.parent_ids = vec![10]; // explicit parent RA

        let mut del_entry = GdEntry::del(60, "chr1", 10, 20);
        del_entry.parent_ids = vec![20]; // explicit parent MC

        let all_records = vec![
            ra_entry.clone(),
            mc_entry.clone(),
            jc_entry.clone(),
            mob_entry.clone(),
            snp_entry.clone(),
            del_entry.clone(),
        ];

        // 1. SNP resolves RA evidence
        let ev_snp = EvidenceReferences::resolve_for_group(&[&snp_entry], Some(&all_records));
        assert_eq!(
            ev_snp.ra,
            Some(vec![GdEvidenceRef {
                gd_id: 10,
                kind: EvidenceKind::Ra
            }])
        );
        assert_eq!(ev_snp.mc, None);
        assert_eq!(ev_snp.jc, None);
        assert_eq!(
            ev_snp.pd_attributes.get("pd_ra_depth"),
            Some(&vec!["42.5".to_string()])
        );
        let sum_snp = EvidenceSummary::from_evidence_references(&ev_snp);
        assert_eq!(sum_snp.ra, Some(true));
        assert_eq!(sum_snp.mc, None);
        assert_eq!(sum_snp.jc, None);
        assert_eq!(sum_snp.supporting_reads, Some(35));
        assert_eq!(sum_snp.coverage, Some(42.5));
        assert_eq!(sum_snp.format_brief(), "RA=1;MC=NA;JC=NA");

        // 2. DEL resolves MC evidence
        let ev_del = EvidenceReferences::resolve_for_group(&[&del_entry], Some(&all_records));
        assert_eq!(ev_del.ra, None);
        assert_eq!(
            ev_del.mc,
            Some(vec![GdEvidenceRef {
                gd_id: 20,
                kind: EvidenceKind::Mc
            }])
        );
        assert_eq!(ev_del.jc, None);
        let sum_del = EvidenceSummary::from_evidence_references(&ev_del);
        assert_eq!(sum_del.ra, None);
        assert_eq!(sum_del.mc, Some(true));
        assert_eq!(sum_del.jc, None);
        assert_eq!(sum_del.format_brief(), "RA=NA;MC=1;JC=NA");

        // 3. JC entry self-references JC evidence
        let ev_jc = EvidenceReferences::resolve_for_group(&[&jc_entry], Some(&all_records));
        assert_eq!(
            ev_jc.jc,
            Some(vec![GdEvidenceRef {
                gd_id: 30,
                kind: EvidenceKind::Jc
            }])
        );
        let sum_jc = EvidenceSummary::from_evidence_references(&ev_jc);
        assert_eq!(sum_jc.jc, Some(true));
        assert_eq!(sum_jc.supporting_reads, Some(15));
        assert_eq!(sum_jc.format_brief(), "RA=NA;MC=NA;JC=1");

        // 4. MOB with JC parent_ids resolves JC evidence
        let ev_mob = EvidenceReferences::resolve_for_group(&[&mob_entry], Some(&all_records));
        assert_eq!(
            ev_mob.jc,
            Some(vec![GdEvidenceRef {
                gd_id: 30,
                kind: EvidenceKind::Jc
            }])
        );
        let sum_mob = EvidenceSummary::from_evidence_references(&ev_mob);
        assert_eq!(sum_mob.jc, Some(true));
        assert_eq!(sum_mob.format_brief(), "RA=NA;MC=NA;JC=1");

        // 5. Multi-evidence event (union across group)
        let ev_all = EvidenceReferences::resolve_for_group(
            &[&snp_entry, &del_entry, &jc_entry],
            Some(&all_records),
        );
        let sum_all = EvidenceSummary::from_evidence_references(&ev_all);
        assert_eq!(sum_all.ra, Some(true));
        assert_eq!(sum_all.mc, Some(true));
        assert_eq!(sum_all.jc, Some(true));
        assert_eq!(sum_all.format_brief(), "RA=1;MC=1;JC=1");
    }

    #[test]
    fn evidence_conservative_failures() {
        use crate::audit::EvidenceSummary;

        // 1. Absent RA: SNP without RA parent or pd_ra_* metrics
        let mut snp_no_ev = GdEntry::snp(1, "chr1", 5, "C");
        snp_no_ev.attrs.insert("frequency".into(), "0.99".into()); // non-pd_* attr
        snp_no_ev.attrs.insert("gene".into(), "test".into());

        let ev_absent = EvidenceReferences::resolve_for_group(&[&snp_no_ev], None);
        assert_eq!(ev_absent.ra, None);
        assert_eq!(ev_absent.mc, None);
        assert_eq!(ev_absent.jc, None);
        assert!(
            ev_absent.pd_attributes.is_empty(),
            "non-pd_* attrs must be excluded"
        );
        let sum_absent = EvidenceSummary::from_evidence_references(&ev_absent);
        assert_eq!(sum_absent.ra, None);
        assert_eq!(sum_absent.mc, None);
        assert_eq!(sum_absent.jc, None);
        assert_eq!(sum_absent.format_brief(), "RA=NA;MC=NA;JC=NA");

        // 2. Dangling parent stays in unresolved_parent_ids
        let mut entry_dangling = GdEntry::snp(2, "chr1", 6, "T");
        entry_dangling.parent_ids = vec![999, 888];
        let ev_dangling = EvidenceReferences::resolve_for_group(&[&entry_dangling], Some(&[]));
        assert_eq!(ev_dangling.unresolved_parent_ids, vec![888, 999]);
        assert_eq!(ev_dangling.ra, None);
        assert_eq!(ev_dangling.mc, None);
        assert_eq!(ev_dangling.jc, None);

        // 3. Parent pointing to a non-evidence record (e.g. another SNP) is unresolved
        let other_snp = GdEntry::snp(50, "chr1", 10, "G");
        let mut entry_with_snp_parent = GdEntry::snp(3, "chr1", 7, "A");
        entry_with_snp_parent.parent_ids = vec![50];
        let ev_non_ev_parent =
            EvidenceReferences::resolve_for_group(&[&entry_with_snp_parent], Some(&[other_snp]));
        assert_eq!(ev_non_ev_parent.unresolved_parent_ids, vec![50]);
        assert_eq!(ev_non_ev_parent.ra, None);

        // 4. Malformed pd_* attributes remain raw strings, but supporting_reads/coverage are None
        let mut entry_malformed = GdEntry::snp(4, "chr1", 8, "C");
        entry_malformed
            .attrs
            .insert("pd_ra_support_reads".into(), "not_an_int".into());
        entry_malformed
            .attrs
            .insert("pd_ra_depth".into(), "not_a_float".into());
        let ev_malformed = EvidenceReferences::resolve_for_group(&[&entry_malformed], None);
        assert_eq!(
            ev_malformed.pd_attributes.get("pd_ra_support_reads"),
            Some(&vec!["not_an_int".to_string()])
        );
        assert_eq!(
            ev_malformed.pd_attributes.get("pd_ra_depth"),
            Some(&vec!["not_a_float".to_string()])
        );
        let sum_malformed = EvidenceSummary::from_evidence_references(&ev_malformed);
        assert_eq!(sum_malformed.supporting_reads, None);
        assert_eq!(sum_malformed.coverage, None);
    }
}
