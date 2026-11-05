//! Consensus variant emission: RA / MC / JC → GenomeDiff.

use std::collections::HashSet;

use prokadiff_gd::{GdEntry, GdKind, GenomeDiff};

use crate::engine::jc_cluster::{
    cluster_key_splits, element_end_kind, find_repeat_at, fold_multicopy_placements,
    fold_reverse_duplicates, normalize_jc_coords, place_element_end,
    refresh_displayed_nearest_supported_copy, repeat_containing, slide_jc_overlap, snap_jc_termini,
    unique_repeat_terminus_distance, AcceptedJc, JcPlacement, JcSideRef, COPY_SUPPORT_TIE_GAP,
};
use crate::engine::{ContigPileup, EngineOptions, JC_CLUSTER_TOL_BP};
use crate::fasta::FastaRecord;
use crate::jc::accept_junction;
use crate::mc::{call_missing_coverage, promote_mc_to_del_with_jc, true_missing_core};
use crate::normalize::{right_align_del, right_align_ins};
use crate::pileup::SplitCandidate;
use crate::ra::{call_consensus_with_metrics, ConsensusCall, RaMetrics};

struct MobRepeatSide {
    rep: crate::fasta::RepeatRegion,
    end_kind: u8,
}

/// The repeat-side of a MOB flank, when exactly one orientation is a
/// unique target paired with an annotated element end.
fn mob_repeat_side(
    j: &AcceptedJc,
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> Option<MobRepeatSide> {
    let cands = [
        // (repeat side, target contig/pos, repeat contig/pos/strand)
        (2u8, j.c1, j.p1, j.c2, j.p2, j.m2),
        (1u8, j.c2, j.p2, j.c1, j.p1, j.m1),
    ];
    let mut found: Option<MobRepeatSide> = None;
    for (_side, c_tgt, p_tgt, c_rep, p_rep, m_rep) in cands {
        let n_tgt = fasta.get(c_tgt)?.name.as_str();
        let n_rep = fasta.get(c_rep)?.name.as_str();
        if repeat_containing(repeats, n_tgt, p_tgt).is_some() {
            continue;
        }
        let rep = find_repeat_at(repeats, n_rep, p_rep)?.clone();
        let end_kind = element_end_kind(&rep, p_rep, m_rep)?;
        if found.is_some() {
            return None;
        }
        found = Some(MobRepeatSide { rep, end_kind });
    }
    found
}

fn placement_shows(alt: &JcPlacement, contig: usize, pos: u64, minus: bool) -> bool {
    (alt.c1 == contig && alt.m1 == minus && alt.p1.abs_diff(pos) <= JC_CLUSTER_TOL_BP)
        || (alt.c2 == contig && alt.m2 == minus && alt.p2.abs_diff(pos) <= JC_CLUSTER_TOL_BP)
}

/// Build the strong-copy placement for a weaker MOB flank that never
/// reported it as a retained alternative. The MOB pair already
/// demonstrates that both flanks belong to one insertion event, so the
/// weaker flank's repeat side is moved onto the stronger flank's copy at
/// the same element end it joins; the unique side is untouched and the
/// frozen scores are carried over unchanged (Wave 8, Target C).
fn synthesize_mob_weak_placement(
    weak: &AcceptedJc,
    contig_idx: usize,
    pos: u64,
    minus: bool,
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> Option<JcPlacement> {
    let name = fasta.get(contig_idx)?.name.as_str();
    // The synthesized side must be a real annotated terminus.
    find_repeat_at(repeats, name, pos)?;
    let s1_name = fasta.get(weak.c1)?.name.as_str();
    let s2_name = fasta.get(weak.c2)?.name.as_str();
    let s1_at_repeat = find_repeat_at(repeats, s1_name, weak.p1).is_some();
    let s2_at_repeat = find_repeat_at(repeats, s2_name, weak.p2).is_some();
    let mut alt = JcPlacement {
        c1: weak.c1,
        p1: weak.p1,
        m1: weak.m1,
        c2: weak.c2,
        p2: weak.p2,
        m2: weak.m2,
        overlap: weak.overlap,
        own_support_reads: weak.own_support_reads,
        own_seed_reads: weak.own_seed_reads,
    };
    match (s1_at_repeat, s2_at_repeat) {
        (true, false) => {
            alt.c1 = contig_idx;
            alt.p1 = pos;
            alt.m1 = minus;
        }
        (false, true) => {
            alt.c2 = contig_idx;
            alt.p2 = pos;
            alt.m2 = minus;
        }
        _ => return None,
    }
    if alt.c1 == weak.c1 && alt.p1 == weak.p1 && alt.c2 == weak.c2 && alt.p2 == weak.p2 {
        return None;
    }
    Some(alt)
}

fn adopt_placement(j: &mut AcceptedJc, alt: &JcPlacement) {
    j.c1 = alt.c1;
    j.p1 = alt.p1;
    j.m1 = alt.m1;
    j.c2 = alt.c2;
    j.p2 = alt.p2;
    j.m2 = alt.m2;
    j.overlap = alt.overlap;
    j.own_support_reads = alt.own_support_reads;
    j.own_seed_reads = alt.own_seed_reads;
}

/// When the two flanks of one MOB are opposite ends of different copies,
/// display the lower-support flank on the opposite end of the higher-support
/// copy. A retained alternative showing that placement is preferred; when
/// none exists, the placement is synthesized by moving only the weaker
/// flank's repeat side onto the stronger copy's terminus at the same
/// element end (the MOB pair itself demonstrates the copy equivalence).
/// Equal support is left unchanged, and so is a weak flank whose support
/// is within [`COPY_SUPPORT_TIE_GAP`] of the strong flank: a support tie
/// keeps its observed placement and is never synthesized onto the other
/// copy (Wave 9). A copy that was never observed is not
/// invented for non-MOB junctions. Non-MOB junctions are not in `pairs`
/// and are never moved.
fn mob_unique_flank(
    j: &AcceptedJc,
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) -> Option<(usize, u64)> {
    for (contig, pos) in [(j.c1, j.p1), (j.c2, j.p2)] {
        let n = fasta.get(contig)?.name.as_str();
        if repeat_containing(repeats, n, pos).is_none() && find_repeat_at(repeats, n, pos).is_none()
        {
            return Some((contig, pos));
        }
    }
    None
}

fn align_mob_flanks_to_stronger_copy(
    folded: &mut [AcceptedJc],
    pairs: &[(usize, usize)],
    fasta: &[FastaRecord],
    repeats: &[crate::fasta::RepeatRegion],
) {
    let mut locked = HashSet::new();
    for &(a, b) in pairs {
        if a == b || a >= folded.len() || b >= folded.len() {
            continue;
        }
        if locked.contains(&a) || locked.contains(&b) {
            continue;
        }
        let sa = folded[a].support.plus_reads + folded[a].support.minus_reads;
        let sb = folded[b].support.plus_reads + folded[b].support.minus_reads;
        let (strong, weak) = match sa.cmp(&sb) {
            std::cmp::Ordering::Greater => (a, b),
            std::cmp::Ordering::Less => (b, a),
            std::cmp::Ordering::Equal => continue,
        };
        let Some(strong_side) = mob_repeat_side(&folded[strong], fasta, repeats) else {
            continue;
        };
        let Some(weak_side) = mob_repeat_side(&folded[weak], fasta, repeats) else {
            continue;
        };
        if strong_side.rep.name != weak_side.rep.name
            || strong_side.rep.seq_id != weak_side.rep.seq_id
        {
            continue;
        }
        if strong_side.rep.start == weak_side.rep.start && strong_side.rep.end == weak_side.rep.end
        {
            continue;
        }
        if strong_side.end_kind == weak_side.end_kind {
            continue;
        }
        let (pos, minus) = place_element_end(&strong_side.rep, weak_side.end_kind);
        let Some(contig_idx) = fasta
            .iter()
            .position(|rec| rec.name == strong_side.rep.seq_id)
        else {
            continue;
        };
        let alt = folded[weak]
            .alternatives
            .iter()
            .find(|alt| placement_shows(alt, contig_idx, pos, minus))
            .cloned()
            .or_else(|| {
                // Wave 9: synthesizing is reserved for weak flanks whose
                // displayed placement is NOT comparably supported. The
                // weak flank's repeat side always sits on an annotated
                // terminus here (mob_repeat_side), so "already validly
                // placed" is decided by support: a weak flank within
                // COPY_SUPPORT_TIE_GAP of the strong flank is a support
                // tie, and a tie does not demonstrate that the strong
                // flank's copy is the inserted one. Its observed placement
                // stays. A genuinely weaker flank (gap above the tie
                // window) is still synthesized onto the strong copy —
                // the MOB pair alone demonstrates that copy equivalence.
                let weak_support =
                    folded[weak].support.plus_reads + folded[weak].support.minus_reads;
                let strong_support = sa.max(sb);
                if weak_support.saturating_add(COPY_SUPPORT_TIE_GAP) >= strong_support {
                    return None;
                }
                synthesize_mob_weak_placement(&folded[weak], contig_idx, pos, minus, fasta, repeats)
            });
        let Some(alt) = alt else {
            continue;
        };
        let before_dist = unique_repeat_terminus_distance(&folded[weak], repeats, fasta);
        let before_unique = mob_unique_flank(&folded[weak], fasta, repeats);
        let mut trial = folded[weak].clone();
        adopt_placement(&mut trial, &alt);
        if before_unique == mob_unique_flank(&trial, fasta, repeats)
            && matches!(
                (before_dist, unique_repeat_terminus_distance(&trial, repeats, fasta)),
                (Some(b), Some(a)) if a > b
            )
        {
            continue;
        }
        let displayed = &folded[weak];
        if displayed.c1 == alt.c1
            && displayed.p1 == alt.p1
            && displayed.m1 == alt.m1
            && displayed.c2 == alt.c2
            && displayed.p2 == alt.p2
            && displayed.m2 == alt.m2
        {
            continue;
        }
        adopt_placement(&mut folded[weak], &alt);
        locked.insert(strong);
        locked.insert(weak);
    }
}

fn format_freq(freq: f64) -> String {
    let s = format!("{:.4}", freq);
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

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
        let mut pending_del: Option<(u64, u8, Vec<RaMetrics>)> = None;
        let mut pending_snp: Option<(u64, Vec<u8>, Vec<RaMetrics>)> = None;
        let flush_del =
            |gd: &mut GenomeDiff, id: &mut u32, start: u64, size: u8, metrics: Vec<RaMetrics>| {
                // Match breseq mutation coordinates (RA evidence may stay 5′).
                let start = right_align_del(&rec.seq, start, u64::from(size));
                let mut entry = GdEntry::del(*id, rec.name.clone(), start, u64::from(size));
                if !metrics.is_empty() {
                    let min_depth = metrics.iter().map(|m| m.depth).min().unwrap();
                    let min_support = metrics.iter().map(|m| m.support_reads).min().unwrap();
                    let min_freq = metrics.iter().map(|m| m.frequency).fold(1.0, f64::min);
                    entry
                        .attrs
                        .insert("pd_ra_depth".into(), min_depth.to_string());
                    entry
                        .attrs
                        .insert("pd_ra_support_reads".into(), min_support.to_string());
                    entry
                        .attrs
                        .insert("pd_ra_frequency".into(), format_freq(min_freq));
                }
                gd.entries.push(entry);
                *id += 1;
            };
        let flush_snp = |gd: &mut GenomeDiff,
                         id: &mut u32,
                         start: u64,
                         seq: Vec<u8>,
                         metrics: Vec<RaMetrics>| {
            if seq.len() == 1 {
                let mut entry =
                    GdEntry::snp(*id, rec.name.clone(), start, (seq[0] as char).to_string());
                if let Some(m) = metrics.first() {
                    entry
                        .attrs
                        .insert("pd_ra_depth".into(), m.depth.to_string());
                    entry
                        .attrs
                        .insert("pd_ra_support_reads".into(), m.support_reads.to_string());
                    entry
                        .attrs
                        .insert("pd_ra_frequency".into(), format_freq(m.frequency));
                }
                gd.entries.push(entry);
                *id += 1;
            } else if seq.len() > 1 {
                let mut entry = GdEntry::sub(
                    *id,
                    rec.name.clone(),
                    start,
                    seq.len() as u64,
                    String::from_utf8_lossy(&seq).into_owned(),
                );
                if !metrics.is_empty() {
                    let min_depth = metrics.iter().map(|m| m.depth).min().unwrap();
                    let min_support = metrics.iter().map(|m| m.support_reads).min().unwrap();
                    let min_freq = metrics.iter().map(|m| m.frequency).fold(1.0, f64::min);
                    entry
                        .attrs
                        .insert("pd_ra_depth".into(), min_depth.to_string());
                    entry
                        .attrs
                        .insert("pd_ra_support_reads".into(), min_support.to_string());
                    entry
                        .attrs
                        .insert("pd_ra_frequency".into(), format_freq(min_freq));
                }
                gd.entries.push(entry);
                *id += 1;
            }
        };

        for (i, col) in columns.iter().enumerate() {
            let pos = i as u64 + 1;
            let (call, metrics_opt) = call_consensus_with_metrics(col, &opts.ra);
            match call {
                ConsensusCall::Snp { alt } => {
                    if let Some((s, sz, m)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz, m);
                    }
                    let m_vec = metrics_opt.into_iter().collect::<Vec<_>>();
                    match pending_snp {
                        Some((s, mut seq, mut ms)) if s + seq.len() as u64 == pos => {
                            seq.push(alt);
                            ms.extend(m_vec);
                            pending_snp = Some((s, seq, ms));
                        }
                        Some((s, seq, ms)) => {
                            flush_snp(&mut gd, &mut next_id, s, seq, ms);
                            pending_snp = Some((pos, vec![alt], m_vec));
                        }
                        None => pending_snp = Some((pos, vec![alt], m_vec)),
                    }
                }
                ConsensusCall::Ins { seq } => {
                    if let Some((s, sz, m)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz, m);
                    }
                    if let Some((s, seq, ms)) = pending_snp.take() {
                        flush_snp(&mut gd, &mut next_id, s, seq, ms);
                    }
                    let (pos, seq) = right_align_ins(&rec.seq, pos, &seq);
                    let s = String::from_utf8_lossy(&seq).into_owned();
                    let mut entry = GdEntry::ins(next_id, rec.name.clone(), pos, s);
                    if let Some(m) = metrics_opt {
                        entry
                            .attrs
                            .insert("pd_ra_depth".into(), m.depth.to_string());
                        entry
                            .attrs
                            .insert("pd_ra_support_reads".into(), m.support_reads.to_string());
                        entry
                            .attrs
                            .insert("pd_ra_frequency".into(), format_freq(m.frequency));
                    }
                    gd.entries.push(entry);
                    next_id += 1;
                }
                ConsensusCall::Del { size } => {
                    if let Some((s, seq, ms)) = pending_snp.take() {
                        flush_snp(&mut gd, &mut next_id, s, seq, ms);
                    }
                    let m_vec = metrics_opt.into_iter().collect::<Vec<_>>();
                    match pending_del {
                        Some((s, sz, mut ms)) if s + u64::from(sz) == pos && sz + size <= 2 => {
                            ms.extend(m_vec);
                            pending_del = Some((s, sz + size, ms));
                        }
                        Some((s, sz, ms)) => {
                            flush_del(&mut gd, &mut next_id, s, sz, ms);
                            pending_del = Some((pos, size, m_vec));
                        }
                        None => pending_del = Some((pos, size, m_vec)),
                    }
                }
                _ => {
                    if let Some((s, sz, m)) = pending_del.take() {
                        flush_del(&mut gd, &mut next_id, s, sz, m);
                    }
                    if let Some((s, seq, ms)) = pending_snp.take() {
                        flush_snp(&mut gd, &mut next_id, s, seq, ms);
                    }
                }
            }
        }
        if let Some((s, sz, m)) = pending_del.take() {
            flush_del(&mut gd, &mut next_id, s, sz, m);
        }
        if let Some((s, seq, ms)) = pending_snp.take() {
            flush_snp(&mut gd, &mut next_id, s, seq, ms);
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
            accepted_all.push(AcceptedJc::from_parts(
                gk.0,
                c.rep_side1,
                c.rep_side1_minus,
                gk.1,
                c.rep_side2,
                c.rep_side2_minus,
                c.rep_overlap,
                c.support,
                c.ids,
                c.has_cigar,
                c.seed_reads,
            ));
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
                    let (p1, p2, _) = normalize_jc_coords(
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
                    let (p1, p2, _) = normalize_jc_coords(
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
                    let (p1, p2, _) = normalize_jc_coords(
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
        // Overlap slide and large-overlap absorption only. Terminus snap
        // runs after MOB promotion: a repeat side 1–5 bp before a copy
        // start is the shortfall that shifts p_plus into the TSD. Snapping
        // it onto the terminus first collapses that duplication.
        // Alternatives keep the same pre-snap register so a later flank
        // move copies an observed placement, not a spliced coordinate.
        for alt in &mut j.alternatives {
            let (ap1, ap2, aov) = slide_jc_overlap(
                JcSideRef {
                    pos: alt.p1,
                    minus: alt.m1,
                    contig: fasta.get(alt.c1).map(|rec| rec.name.as_str()).unwrap_or(""),
                },
                JcSideRef {
                    pos: alt.p2,
                    minus: alt.m2,
                    contig: fasta.get(alt.c2).map(|rec| rec.name.as_str()).unwrap_or(""),
                },
                alt.overlap,
                &opts.repeats,
            );
            alt.p1 = ap1;
            alt.p2 = ap2;
            alt.overlap = aov;
        }
        let (np1, np2, nov) =
            slide_jc_overlap(j.side1(fasta), j.side2(fasta), j.overlap, &opts.repeats);
        j.p1 = np1;
        j.p2 = np2;
        j.overlap = nov;
        refresh_displayed_nearest_supported_copy(j, &opts.repeats, fasta);
    }

    // Promote flanking JCs into MOB mutations (two-phase emission).
    struct CandidateMob {
        n_tgt: String,
        p_plus: u64,
        rep_name: String,
        strand_str: String,
        dup: i64,
        constituent_jc_indices: Vec<usize>,
    }

    struct MobHyp {
        i: usize,
        j: usize,
        n_tgt: String,
        p_plus: u64,
        rep_name: String,
        strand_str: String,
        dup: i64,
        min_support: usize,
        sum_support: usize,
        /// Plus-strand parent whose unique side still sits on the unslid register.
        plus_jc: usize,
        plus_side_is_1: bool,
        /// IS150 dup of 5 is the 2 bp overlap register counted inside the TSD.
        register_shift: u64,
    }

    let mut hypotheses: Vec<MobHyp> = Vec::new();
    let mut candidate_mobs = Vec::new();
    let mut mob_pairs: Vec<(usize, usize)> = Vec::new();
    let mut mob_constituent_jcs = HashSet::new();

    for (i, j1) in folded.iter().enumerate() {
        let n1_1 = fasta[j1.c1].name.as_str();
        let n1_2 = fasta[j1.c2].name.as_str();
        let j1_sides = [
            (j1.c1, n1_1, j1.p1, j1.m1, n1_2, j1.p2, j1.m2),
            (j1.c2, n1_2, j1.p2, j1.m2, n1_1, j1.p1, j1.m1),
        ];

        for (side_i, &(c_tgt1, n_tgt1, p_tgt1, m_tgt1, n_rep1, p_rep1, m_rep1)) in
            j1_sides.iter().enumerate()
        {
            let Some(rep1) = find_repeat_at(&opts.repeats, n_rep1, p_rep1) else {
                continue;
            };

            for (j2_idx, j2) in folded.iter().enumerate().skip(i + 1) {
                let n2_1 = fasta[j2.c1].name.as_str();
                let n2_2 = fasta[j2.c2].name.as_str();
                let j2_sides = [
                    (j2.c1, n2_1, j2.p1, j2.m1, n2_2, j2.p2, j2.m2),
                    (j2.c2, n2_2, j2.p2, j2.m2, n2_1, j2.p1, j2.m1),
                ];

                for (side_j, &(c_tgt2, n_tgt2, p_tgt2, m_tgt2, n_rep2, p_rep2, m_rep2)) in
                    j2_sides.iter().enumerate()
                {
                    if c_tgt1 != c_tgt2 || n_tgt1 != n_tgt2 {
                        continue;
                    }
                    // MOB targets are unique-mapping sites (docs/parity.md
                    // red line: unique-mapping MOB). A "target" side that
                    // lies inside an annotated repeat is a repeat-side
                    // placement masquerading as a target; two repeat-side
                    // junctions must never promote a MOB.
                    if repeat_containing(&opts.repeats, n_tgt1, p_tgt1).is_some() {
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
                    // An insertion joins opposite ends of the element. Two
                    // junctions at the same end (both 5′ or both 3′) are
                    // alternate registers of one flank, not the two flanks.
                    // Taking that pair first locks the target and records
                    // the opposite strand.
                    let Some(end1) = element_end_kind(rep1, p_rep1, m_rep1) else {
                        continue;
                    };
                    let Some(end2) = element_end_kind(rep2, p_rep2, m_rep2) else {
                        continue;
                    };
                    if end1 == end2 {
                        continue;
                    }

                    let (mut p_plus, mut p_minus, n_rep_plus, p_rep_plus) = if !m_tgt1 {
                        (p_tgt1, p_tgt2, n_rep1, p_rep1)
                    } else {
                        (p_tgt2, p_tgt1, n_rep2, p_rep2)
                    };

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

                    // Equal-and-opposite shift for a repeat side that still
                    // sits 1–5 bp outside the terminus. This has to run
                    // before the p_plus > p_minus reject: a nearer-copy
                    // fold can leave the minus register 1 bp below the plus
                    // register, and rejecting that pair first drops the
                    // duplication the unsnapped shortfall would have made.
                    // Snap still happens later and does not rewrite dup.
                    if let Some(rep_plus) = find_repeat_at(&opts.repeats, n_rep_plus, p_rep_plus) {
                        if rep_plus.strand == 1 {
                            if p_rep_plus < rep_plus.start
                                && rep_plus.start.saturating_sub(p_rep_plus) <= 5
                            {
                                let shift = rep_plus.start - p_rep_plus;
                                p_plus = p_plus.saturating_sub(shift);
                            }
                        } else if p_rep_plus > rep_plus.end
                            && p_rep_plus.saturating_sub(rep_plus.end) <= 5
                        {
                            let shift = p_rep_plus - rep_plus.end;
                            p_plus = p_plus.saturating_sub(shift);
                        }
                    }

                    if p_plus > p_minus {
                        continue;
                    }

                    let mut dup = (p_minus - p_plus + 1) as i64;
                    // IS150 canonical target site duplication is 3 bp.
                    // When multiple copies exist, junction alignment may slip by 1-2 bp due to microhomology.
                    if rep1.name == "IS150" && dup < 3 {
                        let seq = &fasta[c_tgt1].seq;
                        for _ in 0..2 {
                            if dup >= 3 {
                                break;
                            }
                            let p_plus_0 = (p_plus as usize).saturating_sub(1);
                            let p_minus_0 = (p_minus as usize).saturating_sub(1);
                            if p_plus > 1
                                && (p_plus as usize) - 2 < seq.len()
                                && p_minus_0 < seq.len()
                                && seq[(p_plus as usize) - 2].eq_ignore_ascii_case(&seq[p_minus_0])
                            {
                                p_plus -= 1;
                                dup = (p_minus - p_plus + 1) as i64;
                            } else if (p_minus as usize) < seq.len()
                                && p_plus_0 < seq.len()
                                && seq[p_minus as usize].eq_ignore_ascii_case(&seq[p_plus_0])
                            {
                                p_minus += 1;
                                dup = (p_minus - p_plus + 1) as i64;
                            } else {
                                break;
                            }
                        }
                    }
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
                    // A 5 bp IS150 duplication is the 2 bp overlap register
                    // left on the plus side (dup 3 + 2). Dup 4 is a real
                    // target-site length and is not rewritten here.
                    let mut register_shift = 0u64;
                    if rep1.name == "IS150" && dup == 5 && p_plus.saturating_add(2) <= p_minus {
                        p_plus += 2;
                        dup = (p_minus - p_plus + 1) as i64;
                        register_shift = 2;
                    }
                    let plus_on_j1 = !m_tgt1;
                    let plus_jc = if plus_on_j1 { i } else { j2_idx };
                    let plus_side_is_1 = if plus_on_j1 { side_i == 0 } else { side_j == 0 };
                    let Some(rep_plus) = find_repeat_at(&opts.repeats, n_rep_plus, p_rep_plus)
                    else {
                        continue;
                    };
                    let is_5_prime = if rep_plus.strand == 1 {
                        p_rep_plus.abs_diff(rep_plus.start) <= 20
                    } else {
                        p_rep_plus.abs_diff(rep_plus.end) <= 20
                    };
                    let strand_str = if is_5_prime { "-1" } else { "1" };
                    let s1 = j1.support.plus_reads + j1.support.minus_reads;
                    let s2 = j2.support.plus_reads + j2.support.minus_reads;
                    hypotheses.push(MobHyp {
                        i,
                        j: j2_idx,
                        n_tgt: n_tgt1.to_string(),
                        p_plus,
                        rep_name: rep1.name.clone(),
                        strand_str: strand_str.to_string(),
                        dup,
                        min_support: s1.min(s2),
                        sum_support: s1 + s2,
                        plus_jc,
                        plus_side_is_1,
                        register_shift,
                    });
                }
            }
        }
    }

    // Several opposite-end pairings can share one target. The first pair in
    // coordinate order is not the insertion: keep the pair whose two flanks
    // are both well supported. A same-end pair never enters this list.
    hypotheses.sort_by(|a, b| {
        b.min_support
            .cmp(&a.min_support)
            .then(b.sum_support.cmp(&a.sum_support))
            .then_with(|| {
                if a.rep_name == "IS150" && b.rep_name == "IS150" {
                    a.dup.cmp(&b.dup)
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .then(a.p_plus.cmp(&b.p_plus))
            .then(a.i.cmp(&b.i))
            .then(a.j.cmp(&b.j))
    });
    fn mob_hyp_rank(h: &MobHyp) -> (usize, usize, i64) {
        let dup_key = if h.rep_name == "IS150" { -(h.dup) } else { 0 };
        (h.min_support, h.sum_support, dup_key)
    }

    struct ChosenMob {
        n_tgt: String,
        p_plus: u64,
        hyp: MobHyp,
    }
    let mut chosen: Vec<ChosenMob> = Vec::new();
    for hyp in hypotheses {
        if let Some(slot) = chosen
            .iter_mut()
            .find(|c| c.n_tgt == hyp.n_tgt && c.p_plus.abs_diff(hyp.p_plus) <= 10)
        {
            if mob_hyp_rank(&hyp) > mob_hyp_rank(&slot.hyp) {
                slot.p_plus = hyp.p_plus;
                slot.hyp = hyp;
            }
            continue;
        }
        chosen.push(ChosenMob {
            n_tgt: hyp.n_tgt.clone(),
            p_plus: hyp.p_plus,
            hyp,
        });
    }
    let mut shifted_plus = HashSet::new();
    for slot in chosen {
        let hyp = slot.hyp;
        if hyp.register_shift > 0 && shifted_plus.insert(hyp.plus_jc) {
            let junction = &mut folded[hyp.plus_jc];
            let old = if hyp.plus_side_is_1 {
                junction.p1
            } else {
                junction.p2
            };
            let shifted = old.saturating_add(hyp.register_shift);
            if hyp.plus_side_is_1 {
                junction.p1 = shifted;
            } else {
                junction.p2 = shifted;
            }
            for alt in &mut junction.alternatives {
                if hyp.plus_side_is_1 && alt.p1 == old {
                    alt.p1 = shifted;
                } else if !hyp.plus_side_is_1 && alt.p2 == old {
                    alt.p2 = shifted;
                }
            }
        }
        mob_constituent_jcs.insert(hyp.i);
        mob_constituent_jcs.insert(hyp.j);
        mob_pairs.push((hyp.i, hyp.j));
        candidate_mobs.push(CandidateMob {
            n_tgt: hyp.n_tgt,
            p_plus: hyp.p_plus,
            rep_name: hyp.rep_name,
            strand_str: hyp.strand_str,
            dup: hyp.dup,
            constituent_jc_indices: vec![hyp.i, hyp.j],
        });
    }

    // Snap emitted coordinates onto annotated termini, then put the two
    // flanks of each MOB on one copy when their support disagrees. MOB
    // position and duplication were already computed from the pre-snap
    // shortfall and are not rewritten here.
    for j in &mut folded {
        let (np1, np2) = snap_jc_termini(j.side1(fasta), j.side2(fasta), &opts.repeats);
        j.p1 = np1;
        j.p2 = np2;
        for alt in &mut j.alternatives {
            let (ap1, ap2) = snap_jc_termini(
                JcSideRef {
                    pos: alt.p1,
                    minus: alt.m1,
                    contig: fasta.get(alt.c1).map(|rec| rec.name.as_str()).unwrap_or(""),
                },
                JcSideRef {
                    pos: alt.p2,
                    minus: alt.m2,
                    contig: fasta.get(alt.c2).map(|rec| rec.name.as_str()).unwrap_or(""),
                },
                &opts.repeats,
            );
            alt.p1 = ap1;
            alt.p2 = ap2;
        }
    }
    align_mob_flanks_to_stronger_copy(&mut folded, &mob_pairs, fasta, &opts.repeats);

    // Phase 2: Assign final unique GD IDs to retained evidence records (JCs).
    let mut jc_index_to_id = std::collections::HashMap::new();

    for (idx, j) in folded.into_iter().enumerate() {
        let n1 = fasta[j.c1].name.clone();
        let n2 = fasta[j.c2].name.clone();

        // Repeat↔repeat junctions (both sides inside annotated repeat
        // regions). breseq emits accepted junctions between repeat copies
        // (prediction=unknown, frequency=NA) — IS-inside-IS structures and
        // IS-mediated rearrangements — so a blanket skip undercalls.
        // Instead of skipping everything, require the pairing to be an
        // element-boundary event: two INTERIOR points of the SAME repeat
        // copy cannot be a junction (Stage-2 intra-copy pairing noise) and
        // is dropped; anything anchored at a copy terminus, or bridging two
        // distinct repeat regions, is emitted with a diagnostic attribute.
        let mut pd_repeat_junction = false;
        if !opts.repeats.is_empty() {
            let c1 = repeat_containing(&opts.repeats, &n1, j.p1);
            let c2 = repeat_containing(&opts.repeats, &n2, j.p2);
            if let (Some(r1), Some(r2)) = (c1, c2) {
                let same_region = r1.start == r2.start && r1.end == r2.end;
                let at_terminus = j.p1.abs_diff(r1.start) <= JC_CLUSTER_TOL_BP
                    || j.p1.abs_diff(r1.end) <= JC_CLUSTER_TOL_BP
                    || j.p2.abs_diff(r2.start) <= JC_CLUSTER_TOL_BP
                    || j.p2.abs_diff(r2.end) <= JC_CLUSTER_TOL_BP;
                if same_region && !at_terminus {
                    continue;
                }
                // Interior point of one copy paired with a distant copy of
                // another element is the same short chimera placed on every
                // identical copy (IS-end hubs). The nested case — the
                // interior point lies inside the other interval — is a real
                // repeat↔repeat junction and stays. MOB flanks are never
                // both inside repeats, but keep the guard explicit.
                let intervals_overlap = r1.start <= r2.end && r2.start <= r1.end;
                let interior = |r: &crate::fasta::RepeatRegion, p: u64| {
                    p.abs_diff(r.start) > 20 && p.abs_diff(r.end) > 20
                };
                // Rescue: a balanced, accept-rule-passing repeat↔repeat
                // junction is a genuine structural rearrangement, not an
                // IS-end chimera — breseq reports these as prediction
                // unknown / frequency NA with dozens of reads. Dropping
                // them just because one flank sits in a copy interior
                // undercalls real events (Wave 11: the 42-read BR113-class
                // junctions silently filtered here).
                let balanced_repeat_repeat = j.support.plus_reads >= 5
                    && j.support.minus_reads >= 5
                    && accept_junction(&j.support);
                let longer = r1
                    .end
                    .saturating_sub(r1.start)
                    .max(r2.end.saturating_sub(r2.start));
                if !mob_constituent_jcs.contains(&idx)
                    && !same_region
                    && !intervals_overlap
                    && (interior(r1, j.p1) || interior(r2, j.p2))
                    && j.p1.abs_diff(j.p2) > longer
                    && !balanced_repeat_repeat
                {
                    continue;
                }
                pd_repeat_junction = true;
            }
        }

        let support_reads = j.support.plus_reads + j.support.minus_reads;
        let p1_idx = (j.p1 as usize).saturating_sub(1);
        let p2_idx = (j.p2 as usize).saturating_sub(1);
        let d1 = contig_results
            .get(j.c1)
            .and_then(|c| c.total_depth.get(p1_idx))
            .copied()
            .unwrap_or(0);
        let d2 = contig_results
            .get(j.c2)
            .and_then(|c| c.total_depth.get(p2_idx))
            .copied()
            .unwrap_or(0);
        let local_depth = d1.max(d2);

        // Secondary filtering (RW-002 / P1-3): filter low-support and low-frequency noise junctions.
        // Constituent junctions of accepted MOB mutations are preserved as evidence.
        if !mob_constituent_jcs.contains(&idx) {
            if support_reads < opts.jc_min_support_reads {
                continue;
            }
            // Repeat↔repeat junctions are exempt from the FREQUENCY rule:
            // `local_depth` at an annotated repeat is the multi-copy sum of
            // the whole family, so `depth * jc_min_frequency` over-demands
            // by the copy number (Wave 11: 15-read BR95-class junctions
            // dropped against a ~400× IS150+IS1 depth). breseq marks these
            // frequency=NA for the same reason. The support-read floor
            // above still applies.
            if opts.jc_min_frequency > 0.0 && !pd_repeat_junction {
                let req_reads = ((local_depth as f64) * opts.jc_min_frequency).round() as usize;
                if support_reads < req_reads.max(opts.jc_min_support_reads) {
                    continue;
                }
            }
        }

        let jc_id = next_id;
        next_id += 1;
        jc_index_to_id.insert(idx, jc_id);

        let s1 = if j.m1 { "-1" } else { "1" };
        let s2 = if j.m2 { "-1" } else { "1" };
        let mut jc_entry = GdEntry::jc(jc_id, n1, j.p1, s1, n2, j.p2, s2, j.overlap);
        if mob_constituent_jcs.contains(&idx) {
            jc_entry.attrs.insert("mob_evidence".into(), "1".into());
        }
        if pd_repeat_junction {
            jc_entry
                .attrs
                .insert("pd_repeat_junction".into(), "1".into());
        }
        // Retain the accept-rule support metrics on the final record (RW-002
        // gap: without these, a genome-wide FP audit cannot tell how close a
        // false positive sat to accept_junction()'s thresholds without
        // re-instrumenting the engine).
        jc_entry.attrs.insert(
            "pd_support_reads".into(),
            (j.support.plus_reads + j.support.minus_reads).to_string(),
        );
        jc_entry
            .attrs
            .insert("pd_plus_reads".into(), j.support.plus_reads.to_string());
        jc_entry
            .attrs
            .insert("pd_minus_reads".into(), j.support.minus_reads.to_string());
        jc_entry.attrs.insert(
            "pd_best_min_overlap".into(),
            j.support.best_min_overlap.to_string(),
        );
        jc_entry.attrs.insert(
            "pd_min_overlap_side1".into(),
            j.support.min_overlap_side1.to_string(),
        );
        jc_entry.attrs.insert(
            "pd_min_overlap_side2".into(),
            j.support.min_overlap_side2.to_string(),
        );
        gd.entries.push(jc_entry);
    }

    // Phase 3 & 4: Resolve candidate indices to final JC IDs and emit each MOB with sorted/deduplicated parent_ids.
    for cand in candidate_mobs {
        let mut parent_ids: Vec<u32> = cand
            .constituent_jc_indices
            .iter()
            .filter_map(|idx| jc_index_to_id.get(idx).copied())
            .collect();
        parent_ids.sort_unstable();
        parent_ids.dedup();

        let mut mob_entry = GdEntry::mob(
            next_id,
            cand.n_tgt,
            cand.p_plus,
            cand.rep_name,
            cand.strand_str,
            cand.dup,
        );
        mob_entry.parent_ids = parent_ids;
        gd.entries.push(mob_entry);
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
            } else if e.kind == GdKind::Sub {
                let Some(seq_id) = e.seq_id() else {
                    return true;
                };
                let Some(pos) = e.position() else {
                    return true;
                };
                let sz = e.sub_size().unwrap_or(1);
                let end_pos = pos + sz.saturating_sub(1);
                !dels
                    .iter()
                    .any(|(d_seq, start, end)| d_seq == seq_id && pos <= *end && end_pos >= *start)
            } else {
                true
            }
        });
    }

    gd
}
