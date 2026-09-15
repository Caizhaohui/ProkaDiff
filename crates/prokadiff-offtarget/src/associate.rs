use crate::model::{MutationOffTargetLink, OffTargetSite};
use prokadiff_gd::{GdEntry, GdKind};

/// Links called mutations to candidate off-target sites based on coordinate distance.
///
/// If distance <= `window_bp`, a `MutationOffTargetLink` record is emitted.
pub fn link_mutations_to_sites(
    mutations: &[GdEntry],
    sites: &[OffTargetSite],
    window_bp: u64,
) -> Vec<MutationOffTargetLink> {
    let mut links = Vec::new();

    for entry in mutations {
        let mut_id = format!("mut_{}", entry.id);
        let mut_kind = entry.kind.as_str().to_string();

        let positions = entry_positions(entry);
        for (m_seq, m_pos) in positions {
            for site in sites {
                if site.seq_id != m_seq {
                    continue;
                }
                let dist = if m_pos < site.start {
                    site.start.saturating_sub(m_pos)
                } else {
                    m_pos.saturating_sub(site.end)
                };

                if dist <= window_bp {
                    links.push(MutationOffTargetLink {
                        mutation_id: mut_id.clone(),
                        site_id: site.site_id.clone(),
                        mutation_type: mut_kind.clone(),
                        mutation_position: m_pos,
                        site_start: site.start,
                        site_end: site.end,
                        distance_to_site: dist,
                        mismatches: site.mismatches,
                        pam: site.pam.clone(),
                        cfd_score: site.cfd_score,
                        association_window: window_bp,
                    });
                }
            }
        }
    }

    links
}

fn entry_positions(entry: &GdEntry) -> Vec<(String, u64)> {
    match entry.kind {
        GdKind::Jc => {
            let mut v = Vec::new();
            if let (Some(s), Some(Ok(p))) = (
                entry.fields.first(),
                entry.fields.get(1).map(|x| x.parse::<u64>()),
            ) {
                v.push((s.clone(), p));
            }
            if let (Some(s), Some(Ok(p))) = (
                entry.fields.get(3),
                entry.fields.get(4).map(|x| x.parse::<u64>()),
            ) {
                v.push((s.clone(), p));
            }
            v
        }
        _ => match (entry.seq_id(), entry.position()) {
            (Some(s), Some(p)) => vec![(s.to_string(), p)],
            _ => Vec::new(),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{BulgeType, Strand};

    #[test]
    fn links_snp_within_window() {
        let entry = GdEntry::snp(1, "chr1", 100, "T");
        let site = OffTargetSite {
            site_id: "site_1".into(),
            seq_id: "chr1".into(),
            start: 120,
            end: 142,
            strand: Strand::Plus,
            guide: "GAGTCCGAGCAGAAGAAGAA".into(),
            target_seq: "GAGTCCGAGCAGAAGAAGAATGG".into(),
            pam: "TGG".into(),
            mismatches: 0,
            bulge_type: BulgeType::None,
            bulge_size: 0,
            search_backend: "rust_exact".into(),
            cfd_score: Some(1.0),
            hsu_score: Some(100.0),
        };

        // Distance from pos 100 to site.start 120 is 20 bp <= 50 bp window
        let links = link_mutations_to_sites(&[entry], &[site], 50);
        assert_eq!(links.len(), 1);
        assert_eq!(links[0].mutation_id, "mut_1");
        assert_eq!(links[0].site_id, "site_1");
        assert_eq!(links[0].distance_to_site, 20);
        assert_eq!(links[0].association_window, 50);
    }

    #[test]
    fn does_not_link_mutation_beyond_window() {
        let entry = GdEntry::snp(1, "chr1", 50, "T");
        let site = OffTargetSite {
            site_id: "site_1".into(),
            seq_id: "chr1".into(),
            start: 120,
            end: 142,
            strand: Strand::Plus,
            guide: "GAGTCCGAGCAGAAGAAGAA".into(),
            target_seq: "GAGTCCGAGCAGAAGAAGAATGG".into(),
            pam: "TGG".into(),
            mismatches: 0,
            bulge_type: BulgeType::None,
            bulge_size: 0,
            search_backend: "rust_exact".into(),
            cfd_score: Some(1.0),
            hsu_score: Some(100.0),
        };

        // Distance is 120 - 50 = 70 > 50
        let links = link_mutations_to_sites(&[entry], &[site], 50);
        assert!(links.is_empty());
    }
}
