//! Multi-clone variant alignment and recurrence classification engine.

use std::collections::{BTreeMap, HashMap};

/// Recurrence status of a post-edit variant across evaluated clones.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RecurrenceStatus {
    /// Variant observed only in a single clone.
    CloneSpecific,
    /// Variant observed in a subset of clones (1 < count < total).
    SharedSubset { count: usize, total: usize },
    /// Variant observed across all evaluated clones.
    SharedAll { total: usize },
}

impl RecurrenceStatus {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::CloneSpecific => "clone_specific",
            Self::SharedSubset { .. } => "shared_subset",
            Self::SharedAll { .. } => "shared_all",
        }
    }
}

/// An observed variant instance from an individual clone.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CloneVariantObservation {
    pub clone_id: String,
    pub event_id: String,
    pub seq_id: String,
    pub position: u64,
    pub end: u64,
    pub gd_type: String,
    pub ref_allele: String,
    pub alt_allele: String,
    pub side2_seq_id: Option<String>,
    pub side2_position: Option<u64>,
    /// Exact canonical wire representation from CanonicalEvent.
    pub canonical_key: Vec<u8>,
}

/// A unified cross-clone variant record in the batch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BatchVariantRecord {
    pub batch_variant_id: String,
    pub seq_id: String,
    pub position: u64,
    pub end: u64,
    pub gd_type: String,
    pub ref_allele: String,
    pub alt_allele: String,
    pub side2_seq_id: Option<String>,
    pub side2_position: Option<u64>,
    pub recurrence: RecurrenceStatus,
    /// Mapping of clone_id -> Option<EventId>
    pub clone_presence: HashMap<String, Option<String>>,
}

/// Group key for identifying equivalent mutation loci across clones based on CanonicalEvent wire bytes.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VariantLocusKey {
    pub canonical_key: Vec<u8>,
}

impl CloneVariantObservation {
    pub fn locus_key(&self) -> VariantLocusKey {
        VariantLocusKey {
            canonical_key: self.canonical_key.clone(),
        }
    }
}

/// Align variants across multiple clones and compute recurrence classifications.
pub fn align_batch_variants(
    clone_ids: &[String],
    observations: &[CloneVariantObservation],
) -> Vec<BatchVariantRecord> {
    let total_clones = clone_ids.len();

    // Map: canonical_key -> (Representative observation metadata, Map of clone_id -> list of event_ids)
    struct GroupVal {
        seq_id: String,
        position: u64,
        end: u64,
        gd_type: String,
        ref_allele: String,
        alt_allele: String,
        side2_seq_id: Option<String>,
        side2_position: Option<u64>,
        presence: BTreeMap<String, Vec<String>>,
    }

    let mut grouped: BTreeMap<Vec<u8>, GroupVal> = BTreeMap::new();

    for obs in observations {
        let entry = grouped
            .entry(obs.canonical_key.clone())
            .or_insert_with(|| GroupVal {
                seq_id: obs.seq_id.clone(),
                position: obs.position,
                end: obs.end,
                gd_type: obs.gd_type.clone(),
                ref_allele: obs.ref_allele.clone(),
                alt_allele: obs.alt_allele.clone(),
                side2_seq_id: obs.side2_seq_id.clone(),
                side2_position: obs.side2_position,
                presence: BTreeMap::new(),
            });

        let events = entry.presence.entry(obs.clone_id.clone()).or_default();
        if !events.contains(&obs.event_id) {
            events.push(obs.event_id.clone());
        }
    }

    let mut records = Vec::with_capacity(grouped.len());
    for (i, (_can_bytes, val)) in grouped.into_iter().enumerate() {
        let count = val.presence.len();
        let recurrence = if total_clones <= 1 || count == 1 {
            RecurrenceStatus::CloneSpecific
        } else if count == total_clones {
            RecurrenceStatus::SharedAll { total: total_clones }
        } else {
            RecurrenceStatus::SharedSubset {
                count,
                total: total_clones,
            }
        };

        let mut clone_presence = HashMap::new();
        for cid in clone_ids {
            let ev = val.presence.get(cid).map(|evs| evs.join(";"));
            clone_presence.insert(cid.clone(), ev);
        }

        let batch_variant_id = format!("BV_{:04}", i + 1);

        records.push(BatchVariantRecord {
            batch_variant_id,
            seq_id: val.seq_id,
            position: val.position,
            end: val.end,
            gd_type: val.gd_type,
            ref_allele: val.ref_allele,
            alt_allele: val.alt_allele,
            side2_seq_id: val.side2_seq_id,
            side2_position: val.side2_position,
            recurrence,
            clone_presence,
        });
    }

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_recurrence_classification_three_clones() {
        let clones = vec!["C01".to_string(), "C02".to_string(), "C03".to_string()];

        let observations = vec![
            // SNP shared across all 3 clones (identical canonical_key)
            CloneVariantObservation {
                clone_id: "C01".to_string(),
                event_id: "DV1_AAA".to_string(),
                seq_id: "chr1".to_string(),
                position: 100,
                end: 100,
                gd_type: "SNP".to_string(),
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                side2_seq_id: None,
                side2_position: None,
                canonical_key: b"SNP:chr1:100:T".to_vec(),
            },
            CloneVariantObservation {
                clone_id: "C02".to_string(),
                event_id: "DV1_BBB".to_string(),
                seq_id: "chr1".to_string(),
                position: 100,
                end: 100,
                gd_type: "SNP".to_string(),
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                side2_seq_id: None,
                side2_position: None,
                canonical_key: b"SNP:chr1:100:T".to_vec(),
            },
            CloneVariantObservation {
                clone_id: "C03".to_string(),
                event_id: "DV1_CCC".to_string(),
                seq_id: "chr1".to_string(),
                position: 100,
                end: 100,
                gd_type: "SNP".to_string(),
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                side2_seq_id: None,
                side2_position: None,
                canonical_key: b"SNP:chr1:100:T".to_vec(),
            },
            // DEL shared between C01 and C02 (SharedSubset)
            CloneVariantObservation {
                clone_id: "C01".to_string(),
                event_id: "DV1_DDD".to_string(),
                seq_id: "chr1".to_string(),
                position: 500,
                end: 505,
                gd_type: "DEL".to_string(),
                ref_allele: "6bp".to_string(),
                alt_allele: "none".to_string(),
                side2_seq_id: None,
                side2_position: None,
                canonical_key: b"DEL:chr1:500:6".to_vec(),
            },
            CloneVariantObservation {
                clone_id: "C02".to_string(),
                event_id: "DV1_EEE".to_string(),
                seq_id: "chr1".to_string(),
                position: 500,
                end: 505,
                gd_type: "DEL".to_string(),
                ref_allele: "6bp".to_string(),
                alt_allele: "none".to_string(),
                side2_seq_id: None,
                side2_position: None,
                canonical_key: b"DEL:chr1:500:6".to_vec(),
            },
            // JC endpoint-swapped representations share canonical_key and must merge to SharedSubset
            CloneVariantObservation {
                clone_id: "C01".to_string(),
                event_id: "DV1_JC1".to_string(),
                seq_id: "chr1".to_string(),
                position: 1000,
                end: 1000,
                gd_type: "JC".to_string(),
                ref_allele: "none".to_string(),
                alt_allele: "chr2:2000:overlap=0".to_string(),
                side2_seq_id: Some("chr2".to_string()),
                side2_position: Some(2000),
                canonical_key: b"JC_CANONICAL_CHR1_1000_CHR2_2000".to_vec(),
            },
            CloneVariantObservation {
                clone_id: "C02".to_string(),
                event_id: "DV1_JC2".to_string(),
                seq_id: "chr1".to_string(),
                position: 1000,
                end: 1000,
                gd_type: "JC".to_string(),
                ref_allele: "none".to_string(),
                alt_allele: "chr2:2000:overlap=0".to_string(),
                side2_seq_id: Some("chr2".to_string()),
                side2_position: Some(2000),
                canonical_key: b"JC_CANONICAL_CHR1_1000_CHR2_2000".to_vec(),
            },
            // CON events at same locus with differing donor regions produce distinct canonical_keys
            CloneVariantObservation {
                clone_id: "C01".to_string(),
                event_id: "DV1_CON1".to_string(),
                seq_id: "chr1".to_string(),
                position: 3000,
                end: 3050,
                gd_type: "CON".to_string(),
                ref_allele: "50bp".to_string(),
                alt_allele: "chr1:4000-4050".to_string(),
                side2_seq_id: None,
                side2_position: None,
                canonical_key: b"CON:chr1:3000:50:chr1:4000-4050".to_vec(),
            },
            CloneVariantObservation {
                clone_id: "C02".to_string(),
                event_id: "DV1_CON2".to_string(),
                seq_id: "chr1".to_string(),
                position: 3000,
                end: 3050,
                gd_type: "CON".to_string(),
                ref_allele: "50bp".to_string(),
                alt_allele: "chr2:5000-5050".to_string(),
                side2_seq_id: None,
                side2_position: None,
                canonical_key: b"CON:chr1:3000:50:chr2:5000-5050".to_vec(),
            },
        ];

        let batch_records = align_batch_variants(&clones, &observations);
        assert_eq!(batch_records.len(), 5);

        // 1. Shared all SNP
        let snp = batch_records.iter().find(|r| r.gd_type == "SNP").unwrap();
        assert_eq!(snp.recurrence, RecurrenceStatus::SharedAll { total: 3 });
        assert_eq!(snp.clone_presence.get("C01").unwrap(), &Some("DV1_AAA".to_string()));
        assert_eq!(snp.clone_presence.get("C02").unwrap(), &Some("DV1_BBB".to_string()));
        assert_eq!(snp.clone_presence.get("C03").unwrap(), &Some("DV1_CCC".to_string()));

        // 2. Shared subset DEL
        let del = batch_records.iter().find(|r| r.gd_type == "DEL").unwrap();
        assert_eq!(
            del.recurrence,
            RecurrenceStatus::SharedSubset { count: 2, total: 3 }
        );

        // 3. Endpoint-swapped JC correctly merged into SharedSubset
        let jc = batch_records.iter().find(|r| r.gd_type == "JC").unwrap();
        assert_eq!(
            jc.recurrence,
            RecurrenceStatus::SharedSubset { count: 2, total: 3 }
        );
        assert_eq!(jc.clone_presence.get("C01").unwrap(), &Some("DV1_JC1".to_string()));
        assert_eq!(jc.clone_presence.get("C02").unwrap(), &Some("DV1_JC2".to_string()));
        assert_eq!(jc.clone_presence.get("C03").unwrap(), &None);

        // 4 & 5. CON events remain separate clone-specific records
        let cons: Vec<_> = batch_records.iter().filter(|r| r.gd_type == "CON").collect();
        assert_eq!(cons.len(), 2);
        for con in cons {
            assert_eq!(con.recurrence, RecurrenceStatus::CloneSpecific);
        }
    }
}
