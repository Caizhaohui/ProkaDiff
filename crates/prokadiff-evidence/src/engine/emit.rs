//! Consensus variant emission: RA / MC / JC → GenomeDiff.

use std::collections::HashSet;

use prokadiff_gd::{GdEntry, GdKind, GenomeDiff};

use crate::engine::jc_cluster::{
    cluster_key_splits, find_repeat_at, fold_multicopy_placements, fold_reverse_duplicates,
    normalize_jc_coords, AcceptedJc,
};
use crate::engine::{ContigPileup, EngineOptions};
use crate::fasta::FastaRecord;
use crate::jc::accept_junction;
use crate::mc::{call_missing_coverage, promote_mc_to_del_with_jc, true_missing_core};
use crate::normalize::{right_align_del, right_align_ins};
use crate::pileup::SplitCandidate;
use crate::ra::{call_consensus, ConsensusCall};

pub(crate) fn emit_from_pileup(
    fasta: &[FastaRecord],
    mut contig_results: Vec<ContigPileup>,
    opts: &EngineOptions,
    extra_splits: &[SplitCandidate],
) -> GenomeDiff {
    for sp in extra_splits {
        if sp.contig_idx < contig_results.len() {
            contig_results[sp.contig_idx].splits.push(sp.clone());
        }
    }

    let mut gd = GenomeDiff::new();
    gd.metadata
        .push(("PROGRAM".into(), "prokadiff-evidence".into()));
    let mut next_id = 1u32;
    // Accepted junctions are collected across ALL contigs and emitted only
    // after redundancy folding: a reverse-direction duplicate keyed on
    // contig B can belong to a physical event first reported on contig A.
    let mut accepted_all: Vec<AcceptedJc> = Vec::new();

    for (idx, rec) in fasta.iter().enumerate() {
        let ContigPileup {
            columns,
            unique_depth,
            total_depth,
            splits,
        } = &contig_results[idx];
        let mut pending_del: Option<(u64, u8)> = None;
        let flush_del = |gd: &mut GenomeDiff, id: &mut u32, start: u64, size: u8| {
            // Match breseq mutation coordinates (RA evidence may stay 5′).
            let start = right_align_del(&rec.seq, start, u64::from(size));
            gd.entries
                .push(GdEntry::del(*id, rec.name.clone(), start, u64::from(size)));
            *id += 1;
        };

        for (i, col) in columns.iter().enumerate() {
            let pos = i as u64 + 1;
            match call_consensus(col, &opts.ra) {
                ConsensusCall::Snp { alt } => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                    gd.entries.push(GdEntry::snp(
                        next_id,
                        rec.name.clone(),
                        pos,
                        (alt as char).to_string(),
                    ));
                    next_id += 1;
                }
                ConsensusCall::Ins { seq } => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                    let (pos, seq) = right_align_ins(&rec.seq, pos, &seq);
                    let s = String::from_utf8_lossy(&seq).into_owned();
                    gd.entries
                        .push(GdEntry::ins(next_id, rec.name.clone(), pos, s));
                    next_id += 1;
                }
                ConsensusCall::Del { size } => match pending_del {
                    Some((s, sz)) if s + u64::from(sz) == pos && sz + size <= 2 => {
                        pending_del = Some((s, sz + size));
                    }
                    Some((s, sz)) => {
                        flush_del(&mut gd, &mut next_id, s, sz);
                        pending_del = Some((pos, size));
                    }
                    None => pending_del = Some((pos, size)),
                },
                _ => {
                    if let Some((s, sz)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz);
                    }
                }
            }
        }
        if let Some((s, sz)) = pending_del.take() {
            flush_del(&mut gd, &mut next_id, s, sz);
        }

        // Junction clustering (replaces exact-key aggregation). Do not key
        // on read strand: plus and minus must aggregate (accept_junction).
        // Group by (contig, side2 contig, overlap sign) — the exact overlap
        // value is register-dependent, but its sign is biologically distinct
        // (deletion-like < 0, blunt = 0, insertion-like > 0). Within a group,
        // sort by (side1, side2) and greedily cluster: a candidate joins the
        // first cluster whose FIRST member is within JC_CLUSTER_TOL_BP on
        // both sides (single-pass against the first member is deterministic;
        // documented in docs/parity.md).
        let mut jc_endpoints_1 = Vec::new();
        for (gk, c) in cluster_key_splits(splits) {
            if !accept_junction(&c.support) {
                continue;
            }
            if gk.0 == idx {
                jc_endpoints_1.push(c.rep_side1);
            }
            if gk.1 == idx {
                jc_endpoints_1.push(c.rep_side2);
            }
            accepted_all.push(AcceptedJc {
                c1: gk.0,
                p1: c.rep_side1,
                m1: c.rep_side1_minus,
                c2: gk.1,
                p2: c.rep_side2,
                m2: c.rep_side2_minus,
                overlap: c.rep_overlap,
                support: c.support,
                ids: c.ids,
                has_cigar: c.has_cigar,
            });
        }

        for span in call_missing_coverage(unique_depth, opts.mc_min_len) {
            if span.len() < 3 {
                continue;
            }
            // Contig ends look like unique-depth 0; breseq reports UN, not DEL.
            if span.start_0 == 0 || span.end_0 == unique_depth.len() {
                continue;
            }
            let start = span.start_1();
            let end = span.end_1();
            let mc_id = next_id;
            next_id += 1;
            gd.entries
                .push(GdEntry::mc(mc_id, rec.name.clone(), start, end, 0, 0));
            if !promote_mc_to_del_with_jc(
                span,
                unique_depth,
                total_depth,
                opts.mc_del_min_len,
                opts.mc_min_len,
                &jc_endpoints_1,
            ) {
                continue;
            }
            let core = true_missing_core(
                span,
                unique_depth,
                total_depth,
                crate::mc::coverage_cutoff(unique_depth),
            )
            .unwrap_or(span);
            let mut del_start = core.start_1();
            let mut del_end = core.end_1();
            let mut mediated_repeat: Option<String> = None;

            // If an accepted JC flanks this missing coverage span,
            // the JC breakpoint defines the exact deletion boundary:
            for j in &accepted_all {
                if j.c1 == idx && j.c2 == idx {
                    let (p1, p2) = normalize_jc_coords(
                        j.side1(fasta),
                        j.side2(fasta),
                        j.overlap,
                        &opts.repeats,
                    );
                    let (norm_jp1, norm_jp2) = if p1 <= p2 { (p1, p2) } else { (p2, p1) };
                    if norm_jp1.abs_diff(span.start_1().saturating_sub(1)) <= 50
                        && norm_jp2.abs_diff(span.end_1().saturating_add(1)) <= 50
                    {
                        del_start = norm_jp1 + 1;
                        del_end = norm_jp2.saturating_sub(1);
                        break;
                    }
                }
            }

            // Also check for repeat-mediated deletions where flanks connect to repeat elements
            let left_anchor = span.start_1().saturating_sub(1);
            let right_anchor = span.end_1().saturating_add(1);
            let left_repeat = opts.repeats.iter().find(|r| {
                r.seq_id == rec.name
                    && (r.end.abs_diff(left_anchor) <= 20 || r.start.abs_diff(left_anchor) <= 20)
            });
            let right_repeat = opts.repeats.iter().find(|r| {
                r.seq_id == rec.name
                    && (r.start.abs_diff(right_anchor) <= 20 || r.end.abs_diff(right_anchor) <= 20)
            });
            if let Some(r) = left_repeat {
                del_start = r.end + 1;
                mediated_repeat = Some(r.name.clone());
                for j in &accepted_all {
                    let (p1, p2) = normalize_jc_coords(
                        j.side1(fasta),
                        j.side2(fasta),
                        j.overlap,
                        &opts.repeats,
                    );
                    if j.c1 == idx
                        && p1.abs_diff(right_anchor) <= 50
                        && p1 > span.end_1().saturating_sub(5)
                    {
                        del_end = p1.saturating_sub(1).min(right_anchor);
                    } else if j.c2 == idx
                        && p2.abs_diff(right_anchor) <= 50
                        && p2 > span.end_1().saturating_sub(5)
                    {
                        del_end = p2.saturating_sub(1).min(right_anchor);
                    }
                }
            } else if let Some(r) = right_repeat {
                del_end = r.start.saturating_sub(1);
                if mediated_repeat.is_none() {
                    mediated_repeat = Some(r.name.clone());
                }
                for j in &accepted_all {
                    let (p1, p2) = normalize_jc_coords(
                        j.side1(fasta),
                        j.side2(fasta),
                        j.overlap,
                        &opts.repeats,
                    );
                    if j.c1 == idx
                        && p1.abs_diff(left_anchor) <= 50
                        && p1 < span.start_1().saturating_add(5)
                    {
                        del_start = (p1 + 1).max(left_anchor);
                    } else if j.c2 == idx
                        && p2.abs_diff(left_anchor) <= 50
                        && p2 < span.start_1().saturating_add(5)
                    {
                        del_start = (p2 + 1).max(left_anchor);
                    }
                }
            }

            let size = del_end.saturating_sub(del_start) + 1;
            let del_start = if mediated_repeat.is_some() {
                del_start
            } else {
                right_align_del(&rec.seq, del_start, size)
            };
            let mut del_entry = GdEntry::del(next_id, rec.name.clone(), del_start, size);
            del_entry.parent_ids = vec![mc_id];
            if let Some(rep_name) = mediated_repeat {
                del_entry.attrs.insert("mediated".into(), rep_name);
            }
            gd.entries.push(del_entry);
            next_id += 1;
        }
    }

    // Redundancy folding (synth_is_mob 2372015: 6 JC emitted for 2
    // physical junctions — 2 multi-copy placements + 2 reverse-direction
    // duplicates). Reverse first: folding copy placements before direction
    // duplicates would strand the copy-side reports. Then emit one JC per
    // physical event, in deterministic coordinate order.
    let mut folded =
        fold_multicopy_placements(fold_reverse_duplicates(accepted_all), &opts.repeats, fasta);
    folded.sort_by_key(|j| (j.c1, j.p1, j.m1, j.c2, j.p2, j.m2));
    for j in &mut folded {
        let (np1, np2) =
            normalize_jc_coords(j.side1(fasta), j.side2(fasta), j.overlap, &opts.repeats);
        j.p1 = np1;
        j.p2 = np2;
        if j.overlap > 0 && j.overlap <= 4 {
            j.overlap = 0;
        }
    }

    // Promote flanking JCs into MOB mutations
    let mut mob_entries = Vec::new();
    let mut used_mobs = HashSet::new();

    for (i, j1) in folded.iter().enumerate() {
        let n1_1 = fasta[j1.c1].name.as_str();
        let n1_2 = fasta[j1.c2].name.as_str();
        let j1_sides = [
            (j1.c1, n1_1, j1.p1, j1.m1, n1_2, j1.p2, j1.m2),
            (j1.c2, n1_2, j1.p2, j1.m2, n1_1, j1.p1, j1.m1),
        ];

        for &(c_tgt1, n_tgt1, p_tgt1, m_tgt1, n_rep1, p_rep1, _m_rep1) in &j1_sides {
            let Some(rep1) = find_repeat_at(&opts.repeats, n_rep1, p_rep1) else {
                continue;
            };

            for j2 in folded.iter().skip(i + 1) {
                let n2_1 = fasta[j2.c1].name.as_str();
                let n2_2 = fasta[j2.c2].name.as_str();
                let j2_sides = [
                    (j2.c1, n2_1, j2.p1, j2.m1, n2_2, j2.p2, j2.m2),
                    (j2.c2, n2_2, j2.p2, j2.m2, n2_1, j2.p1, j2.m1),
                ];

                for &(c_tgt2, n_tgt2, p_tgt2, m_tgt2, n_rep2, p_rep2, _m_rep2) in &j2_sides {
                    if c_tgt1 != c_tgt2 || n_tgt1 != n_tgt2 {
                        continue;
                    }
                    if m_tgt1 == m_tgt2 {
                        continue;
                    }
                    if p_tgt1.abs_diff(p_tgt2) > 10 {
                        continue;
                    }
                    let Some(rep2) = find_repeat_at(&opts.repeats, n_rep2, p_rep2) else {
                        continue;
                    };
                    if rep1.name != rep2.name {
                        continue;
                    }

                    let (p_plus, p_minus, p_rep_plus) = if !m_tgt1 {
                        (p_tgt1, p_tgt2, p_rep1)
                    } else {
                        (p_tgt2, p_tgt1, p_rep2)
                    };

                    if p_plus > p_minus {
                        continue;
                    }
                    let unique_cov_ok = {
                        let u_depth = &contig_results[c_tgt1].unique_depth;
                        let p_idx = (p_plus as usize).saturating_sub(1);
                        let n = u_depth.len();
                        (p_idx.saturating_sub(50)..p_idx.min(n)).any(|pos| u_depth[pos] >= 5)
                            && (p_idx..p_idx.saturating_add(50).min(n)).any(|pos| u_depth[pos] >= 5)
                    };
                    if !unique_cov_ok {
                        continue;
                    }

                    let mut dup = (p_minus - p_plus + 1) as i64;
                    // Biological constraint: IS elements have minimum target duplications (>= 3 bp).
                    // Sub-3 bp duplications are alignment artifacts / unresolved repeat boundaries.
                    if dup < 3 {
                        continue;
                    }
                    // IS150 canonical target site duplication is 3 bp.
                    // If dup == 4 and the endpoints of the target duplication share a 1-bp microhomology
                    // (same base at p_plus and p_minus), the junction register over-counted the TSD by 1 bp.
                    if rep1.name == "IS150" && dup == 4 {
                        let p_plus_0 = (p_plus as usize).saturating_sub(1);
                        let p_minus_0 = (p_minus as usize).saturating_sub(1);
                        let seq = &fasta[c_tgt1].seq;
                        if p_plus_0 < seq.len()
                            && p_minus_0 < seq.len()
                            && seq[p_plus_0].eq_ignore_ascii_case(&seq[p_minus_0])
                        {
                            dup = 3;
                        }
                    }
                    let Some(rep_plus) = find_repeat_at(&opts.repeats, n_tgt1, p_rep_plus) else {
                        continue;
                    };
                    let is_5_prime = if rep_plus.strand == 1 {
                        p_rep_plus.abs_diff(rep_plus.start) <= 20
                    } else {
                        p_rep_plus.abs_diff(rep_plus.end) <= 20
                    };
                    let strand_str = if is_5_prime { "-1" } else { "1" };
                    let key = (n_tgt1.to_string(), p_plus);
                    if used_mobs.insert(key) {
                        mob_entries.push(GdEntry::mob(
                            next_id,
                            n_tgt1,
                            p_plus,
                            rep1.name.as_str(),
                            strand_str,
                            dup,
                        ));
                        next_id += 1;
                    }
                }
            }
        }
    }
    gd.entries.extend(mob_entries);

    for j in folded {
        let s1 = if j.m1 { "-1" } else { "1" };
        let s2 = if j.m2 { "-1" } else { "1" };
        let n1 = fasta[j.c1].name.clone();
        let n2 = fasta[j.c2].name.clone();
        gd.entries
            .push(GdEntry::jc(next_id, n1, j.p1, s1, n2, j.p2, s2, j.overlap));
        next_id += 1;
    }

    // Subsume/mask any SNP that falls inside a called DEL
    let dels: Vec<(String, u64, u64)> = gd
        .entries
        .iter()
        .filter(|e| e.kind == GdKind::Del)
        .filter_map(|e| {
            let seq_id = e.seq_id()?.to_string();
            let pos = e.position()?;
            let sz = e.del_size()?;
            Some((seq_id, pos, pos + sz.saturating_sub(1)))
        })
        .collect();
    if !dels.is_empty() {
        gd.entries.retain(|e| {
            if e.kind == GdKind::Snp {
                let Some(seq_id) = e.seq_id() else {
                    return true;
                };
                let Some(pos) = e.position() else {
                    return true;
                };
                !dels
                    .iter()
                    .any(|(d_seq, start, end)| d_seq == seq_id && pos >= *start && pos <= *end)
            } else {
                true
            }
        });
    }

    gd
}
