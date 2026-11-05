//! Subclonal variant representation and classification aggregation.
//!
//! Milestone M11: Subclonal / Pool-seq Variant Audit Mode.

pub use prokadiff_evidence::subclonal::{
    bh_adjust, binom_upper_tail, calculate_delta_ci, calculate_wilson_ci,
    classify_subclonal_candidate, fisher_two_sided, AlleleCounts, StrandCounts, SubclonalClass,
    SubclonalParams, WilsonCi,
};

/// Authorized small variant kind strictly enforced for subclonal audit mode.
/// Excludes MOB, JC, AMP, and any indels over 2 bp.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AuthorizedSubclonalKind {
    Snp,
    Sub { len: usize },
    Ins { len: usize },
    Del { len: usize },
}

impl AuthorizedSubclonalKind {
    pub fn from_str(kind_str: &str, ref_seq: &str, alt_seq: &str) -> Option<Self> {
        match kind_str.to_ascii_uppercase().as_str() {
            "SNP" => {
                if ref_seq.len() == 1 && alt_seq.len() == 1 && ref_seq != alt_seq {
                    Some(Self::Snp)
                } else {
                    None
                }
            }
            "SUB" => {
                let len = ref_seq.len();
                if (1..=2).contains(&len) && len == alt_seq.len() && ref_seq != alt_seq {
                    Some(Self::Sub { len })
                } else {
                    None
                }
            }
            "INS" => {
                let len = alt_seq.len();
                if (1..=2).contains(&len) {
                    Some(Self::Ins { len })
                } else {
                    None
                }
            }
            "DEL" => {
                let len = ref_seq.len();
                if (1..=2).contains(&len) {
                    Some(Self::Del { len })
                } else {
                    None
                }
            }
            _ => None, // Strictly rejects MOB, JC, AMP, etc.
        }
    }
}

/// Observed candidate small variant prior to statistical evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct RawSubclonalObservation {
    pub contig: String,
    pub pos: u64,
    pub ref_seq: String,
    pub alt_seq: String,
    pub variant_type: String, // "SNP", "SUB", "INS", "DEL"
    pub starter_counts: AlleleCounts,
    pub edited_counts: AlleleCounts,
    pub starter_error_rate: f64,
    pub edited_error_rate: f64,
}

/// Fully evaluated and classified subclonal variant record.
#[derive(Clone, Debug, PartialEq)]
pub struct SubclonalVariant {
    pub contig: String,
    pub pos: u64,
    pub ref_seq: String,
    pub alt_seq: String,
    pub variant_type: String,

    pub starter_counts: AlleleCounts,
    pub edited_counts: AlleleCounts,

    pub starter_af: Option<f64>,
    pub edited_af: Option<f64>,

    pub starter_ci: Option<WilsonCi>,
    pub edited_ci: Option<WilsonCi>,

    pub starter_strand_bias_p: f64,
    pub edited_strand_bias_p: f64,

    pub starter_error_rate: f64,
    pub edited_error_rate: f64,

    pub starter_p_err: f64,
    pub starter_q_err: f64,

    pub edited_p_err: f64,
    pub edited_q_err: f64,

    pub delta_af: f64,
    pub delta_ci_low: f64,
    pub delta_ci_high: f64,
    pub p_diff: f64,
    pub q_diff: f64,

    /// Starter exclusion upper bound at observed depth when starter ALT count is 0.
    pub starter_exclusion_bound: Option<f64>,

    pub classification: SubclonalClass,
    pub flags: Vec<String>,
}

/// Evaluates a collection of candidate small variant observations in batch,
/// performing multiple hypothesis FDR correction (Benjamini-Hochberg),
/// background error floor enforcement, and deterministic sorting by `(contig, pos, ref, alt)`.
/// Excludes unauthorized variant types (MOB, JC, AMP, indels > 2bp), NOT_SUBCLONAL,
/// and subtracted PRE_EXISTING_IN_STARTER from published findings.
pub fn evaluate_subclonal_candidates(
    raw_candidates: Vec<RawSubclonalObservation>,
    params: &SubclonalParams,
) -> Vec<SubclonalVariant> {
    if raw_candidates.is_empty() {
        return Vec::new();
    }

    // Filter to authorized small variant kinds (strictly rejects MOB, JC, AMP, and indels > 2bp)
    let mut valid_candidates: Vec<RawSubclonalObservation> = raw_candidates
        .into_iter()
        .filter(|obs| {
            AuthorizedSubclonalKind::from_str(&obs.variant_type, &obs.ref_seq, &obs.alt_seq).is_some()
        })
        .collect();

    if valid_candidates.is_empty() {
        return Vec::new();
    }

    // Sort deterministically
    valid_candidates.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then_with(|| a.pos.cmp(&b.pos))
            .then_with(|| a.ref_seq.cmp(&b.ref_seq))
            .then_with(|| a.alt_seq.cmp(&b.alt_seq))
    });

    // 1. Calculate per-locus raw p-values for error testing and differential Fisher
    let mut p_err_starter = Vec::with_capacity(valid_candidates.len());
    let mut p_err_edited = Vec::with_capacity(valid_candidates.len());
    let mut p_diff_vec = Vec::with_capacity(valid_candidates.len());

    for obs in &valid_candidates {
        let auth_kind = AuthorizedSubclonalKind::from_str(&obs.variant_type, &obs.ref_seq, &obs.alt_seq)
            .expect("already filtered to authorized kind");
        let floor = match auth_kind {
            AuthorizedSubclonalKind::Snp | AuthorizedSubclonalKind::Sub { .. } => {
                params.error_floor_snp
            }
            AuthorizedSubclonalKind::Ins { .. } | AuthorizedSubclonalKind::Del { .. } => {
                params.error_floor_indel
            }
        };
        let effective_e0_s = obs.starter_error_rate.max(floor);
        let effective_e0_e = obs.edited_error_rate.max(floor);

        let n_s = obs.starter_counts.total_depth();
        let a_s = obs.starter_counts.alt_counts.total();
        let p_s = binom_upper_tail(a_s, n_s, effective_e0_s).unwrap_or(1.0);
        p_err_starter.push(p_s);

        let n_e = obs.edited_counts.total_depth();
        let a_e = obs.edited_counts.alt_counts.total();
        let p_e = binom_upper_tail(a_e, n_e, effective_e0_e).unwrap_or(1.0);
        p_err_edited.push(p_e);

        // Differential 2x2 table: [[A_e, R_e], [A_s, R_s]]
        let r_e = obs.edited_counts.ref_counts.total();
        let r_s = obs.starter_counts.ref_counts.total();
        let diff_table = [[a_e, r_e], [a_s, r_s]];
        let p_diff = fisher_two_sided(diff_table);
        p_diff_vec.push(p_diff);
    }

    // 2. Benjamini-Hochberg FDR adjustment
    let q_err_starter = bh_adjust(&p_err_starter);
    let q_err_edited = bh_adjust(&p_err_edited);
    let q_diff_vec = bh_adjust(&p_diff_vec);

    // 3. Construct SubclonalVariant records
    let mut results = Vec::with_capacity(valid_candidates.len());

    for (i, obs) in valid_candidates.into_iter().enumerate() {
        let n_s = obs.starter_counts.total_depth();
        let a_s = obs.starter_counts.alt_counts.total();
        let af_s = obs.starter_counts.alt_fraction();
        let ci_s = if n_s > 0 {
            calculate_wilson_ci(a_s, n_s).ok()
        } else {
            None
        };

        let n_e = obs.edited_counts.total_depth();
        let a_e = obs.edited_counts.alt_counts.total();
        let af_e = obs.edited_counts.alt_fraction();
        let ci_e = if n_e > 0 {
            calculate_wilson_ci(a_e, n_e).ok()
        } else {
            None
        };

        // Strand bias Fisher exact tests: [[R_f, R_r], [A_f, A_r]]
        let sb_table_s = [
            [
                obs.starter_counts.ref_counts.fwd as u64,
                obs.starter_counts.ref_counts.rev as u64,
            ],
            [
                obs.starter_counts.alt_counts.fwd as u64,
                obs.starter_counts.alt_counts.rev as u64,
            ],
        ];
        let sb_p_s = fisher_two_sided(sb_table_s);

        let sb_table_e = [
            [
                obs.edited_counts.ref_counts.fwd as u64,
                obs.edited_counts.ref_counts.rev as u64,
            ],
            [
                obs.edited_counts.alt_counts.fwd as u64,
                obs.edited_counts.alt_counts.rev as u64,
            ],
        ];
        let sb_p_e = fisher_two_sided(sb_table_e);

        let p_diff = p_diff_vec[i];
        let q_diff = q_diff_vec[i];
        let q_s = q_err_starter[i];
        let q_e = q_err_edited[i];

        let (classification, delta_metrics) = classify_subclonal_candidate(
            &obs.starter_counts,
            &obs.edited_counts,
            q_s,
            q_e,
            p_diff,
            q_diff,
            params,
        );

        let (delta_af, delta_ci_low, delta_ci_high) = match delta_metrics {
            Some(d) => d,
            None => {
                let default_ci = WilsonCi { low: 0.0, high: 0.0 };
                calculate_delta_ci(
                    af_e.unwrap_or(0.0),
                    ci_e.as_ref().unwrap_or(&default_ci),
                    af_s.unwrap_or(0.0),
                    ci_s.as_ref().unwrap_or(&default_ci),
                )
            }
        };

        let starter_exclusion_bound = if a_s == 0 {
            ci_s.map(|ci| ci.high)
        } else {
            None
        };

        // Collect quality flags
        let mut flags = Vec::new();
        if sb_p_e < params.max_strand_bias_pvalue {
            flags.push("STRAND_BIAS".to_string());
        }
        if a_e >= 8
            && (obs.edited_counts.alt_counts.fwd == 0 || obs.edited_counts.alt_counts.rev == 0)
        {
            flags.push("SINGLE_STRAND_ALT".to_string());
        }
        if n_e > 0 && (obs.edited_counts.other as f64 / n_e as f64) > 0.10 {
            flags.push("HIGH_OTHER_FRACTION".to_string());
        }
        if n_e < params.min_depth {
            flags.push("LOW_DEPTH".to_string());
        }
        if let Some(freq) = af_e {
            if (0.70..0.75).contains(&freq) {
                flags.push("NEAR_CONSENSUS".to_string());
            }
        }
        flags.sort();

        results.push(SubclonalVariant {
            contig: obs.contig,
            pos: obs.pos,
            ref_seq: obs.ref_seq,
            alt_seq: obs.alt_seq,
            variant_type: obs.variant_type,
            starter_counts: obs.starter_counts,
            edited_counts: obs.edited_counts,
            starter_af: af_s,
            edited_af: af_e,
            starter_ci: ci_s,
            edited_ci: ci_e,
            starter_strand_bias_p: sb_p_s,
            edited_strand_bias_p: sb_p_e,
            starter_error_rate: obs.starter_error_rate,
            edited_error_rate: obs.edited_error_rate,
            starter_p_err: p_err_starter[i],
            starter_q_err: q_s,
            edited_p_err: p_err_edited[i],
            edited_q_err: q_e,
            delta_af,
            delta_ci_low,
            delta_ci_high,
            p_diff,
            q_diff,
            starter_exclusion_bound,
            classification,
            flags,
        });
    }

    // Exclude NOT_SUBCLONAL and subtracted PRE_EXISTING_IN_STARTER from reportable findings
    results.retain(|r| r.classification.is_reportable());

    results
}

fn calc_window_error(cols: &[prokadiff_evidence::PileupColumn], target_idx: usize) -> f64 {
    let win_start = target_idx.saturating_sub(50);
    let win_end = (target_idx + 51).min(cols.len());
    let mut mismatches = 0u64;
    let mut total_bases = 0u64;

    for (k, col) in cols.iter().enumerate().take(win_end).skip(win_start) {
        if k == target_idx {
            continue;
        }
        for obs in &col.observations {
            if obs.base != b'-' {
                total_bases += 1;
                if !obs.base.eq_ignore_ascii_case(&col.ref_base) {
                    mismatches += 1;
                }
            }
        }
    }

    if total_bases >= 200 {
        (mismatches as f64) / (total_bases as f64)
    } else {
        0.005
    }
}

/// Per-read snapshot of one pileup column.
///
/// M11 Round-3 repair (P1-1 / P1-2): allele eligibility is assigned per read
/// BEFORE aggregation. The previous implementation subtracted aggregate
/// insertion counts from aggregate REF base counts with `saturating_sub`,
/// which let an insertion whose anchor base was repeat-masked (or mismatched)
/// steal unrelated REF observations. `read_id` threading in the pileup
/// builder makes per-read bookkeeping possible instead.
///
/// Invariants relied on here:
/// * one aligned read contributes at most ONE base observation per column
///   (CIGAR walks monotonically), identified by `read_id`;
/// * a read with an L-bp CIGAR deletion yields L consecutive `b'-'`
///   observations carrying the SAME `read_id` at consecutive columns;
/// * insertion entries exist unconditionally at their right-normalized
///   anchor even when the anchor base observation was repeat-masked.
struct ColumnReadIndex {
    /// read_id -> observed base at this column (uppercase A/C/G/T, `b'-'`, or N).
    base_by_read: std::collections::HashMap<u64, u8>,
    /// read_id -> strand of the base observation.
    strand_by_read: std::collections::HashMap<u64, prokadiff_evidence::Strand>,
    /// read_id -> (insertion oligo uppercase, entry strand) pairs at this column.
    ins_by_read: std::collections::HashMap<u64, Vec<(Vec<u8>, prokadiff_evidence::Strand)>>,
}

impl ColumnReadIndex {
    fn new(col: &prokadiff_evidence::PileupColumn) -> Self {
        let mut base_by_read = std::collections::HashMap::new();
        let mut strand_by_read = std::collections::HashMap::new();
        for obs in &col.observations {
            // First observation wins; duplicates for one read_id in one
            // column cannot arise from the builder, this is defensive only.
            base_by_read
                .entry(obs.read_id)
                .or_insert(obs.base.to_ascii_uppercase());
            strand_by_read.entry(obs.read_id).or_insert(obs.strand);
        }
        let mut ins_by_read: std::collections::HashMap<u64, Vec<(Vec<u8>, prokadiff_evidence::Strand)>> =
            std::collections::HashMap::new();
        for ins in &col.insertions {
            ins_by_read
                .entry(ins.read_id)
                .or_default()
                .push((ins.seq.to_ascii_uppercase(), ins.strand));
        }
        Self {
            base_by_read,
            strand_by_read,
            ins_by_read,
        }
    }

    /// Total depth: distinct reads with either a base observation or an
    /// insertion entry at this column. Reads that have both count exactly
    /// once (a repeat-masked anchor still counts via its insertion entry).
    ///
    /// Depth formula (M11 Round-3, P1-2):
    /// `depth(x) = |{read_id : base obs at x} ∪ {read_id : insertion entry at x}|`
    /// `REF + ALT + OTHER == depth(x)` by construction of the per-read split.
    /// Production code derives the same quantity as
    /// `AlleleCounts::total_depth()`; this direct union count backs the
    /// invariant unit tests.
    #[cfg(test)]
    fn depth_union(&self) -> u64 {
        let mut total = self.base_by_read.len() as u64;
        for id in self.ins_by_read.keys() {
            if !self.base_by_read.contains_key(id) {
                total += 1;
            }
        }
        total
    }
}

/// Per-read SNP allele classification at one column.
///
/// * ALT   = reads with base == `alt_b` AND no insertion entry at this column.
/// * REF   = reads with base == ref AND no insertion entry at this column.
/// * OTHER = everything else: mismatches, deletion observations, N bases, and
///   every insertion-carrying read (an insertion-carrying read is not a clean
///   SNP carrier on either allele).
/// * Depth = |REF| + |ALT| + OTHER = `ColumnReadIndex::depth_union`.
fn extract_snp_counts_from_index(idx: &ColumnReadIndex, ref_b: u8, alt_b: u8) -> AlleleCounts {
    let tgt_b = alt_b.to_ascii_uppercase();
    let mut ref_fwd = 0u32;
    let mut ref_rev = 0u32;
    let mut alt_fwd = 0u32;
    let mut alt_rev = 0u32;
    let mut other = 0u32;

    let mut visited: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut ids = idx.base_by_read.keys().copied().chain(idx.ins_by_read.keys().copied());
    for id in std::iter::from_fn(move || ids.next()) {
        if !visited.insert(id) {
            continue;
        }
        let has_ins = idx.ins_by_read.contains_key(&id);
        let base = idx.base_by_read.get(&id).copied();
        if !has_ins && base == Some(tgt_b) {
            match idx.strand_by_read.get(&id) {
                Some(prokadiff_evidence::Strand::Plus) => alt_fwd += 1,
                _ => alt_rev += 1,
            }
        } else if !has_ins && base == Some(ref_b) {
            match idx.strand_by_read.get(&id) {
                Some(prokadiff_evidence::Strand::Plus) => ref_fwd += 1,
                _ => ref_rev += 1,
            }
        } else {
            other += 1;
        }
    }

    AlleleCounts::new(
        StrandCounts::new(ref_fwd, ref_rev),
        StrandCounts::new(alt_fwd, alt_rev),
        other,
    )
}

/// Per-read INS allele classification for candidate oligo `ins_seq` at one column.
///
/// * ALT   = reads with an insertion entry here whose oligo equals `ins_seq`
///   (case-insensitive) — regardless of whether their anchor base observation
///   survived repeat masking; the insertion entry exists unconditionally.
/// * REF   = reads with base == ref AND no insertion entry here (strict set
///   difference — no aggregate subtraction, so pure-REF reads can never be
///   stolen by insertion-carrying reads).
/// * OTHER = reads carrying a different insertion oligo + base mismatches /
///   deletion observations / N bases.
/// * Depth = |{base obs} ∪ {insertion entries}| per read (see depth_union).
fn extract_ins_counts(idx: &ColumnReadIndex, ref_b: u8, ins_seq: &[u8]) -> AlleleCounts {
    let target = ins_seq.to_ascii_uppercase();
    let mut ref_fwd = 0u32;
    let mut ref_rev = 0u32;
    let mut alt_fwd = 0u32;
    let mut alt_rev = 0u32;
    let mut other = 0u32;

    let mut visited: std::collections::HashSet<u64> = std::collections::HashSet::new();

    // Insertion-carrying reads first: they are ALT or OTHER for this candidate.
    for (id, entries) in &idx.ins_by_read {
        visited.insert(*id);
        match entries.iter().find(|(seq, _)| *seq == target) {
            Some((_, strand)) => match strand {
                prokadiff_evidence::Strand::Plus => alt_fwd += 1,
                prokadiff_evidence::Strand::Minus => alt_rev += 1,
            },
            None => other += 1,
        }
    }

    // Reads with only a base observation.
    for id in idx.base_by_read.keys() {
        if !visited.insert(*id) {
            continue;
        }
        let base = idx.base_by_read.get(id).copied();
        if base == Some(ref_b) {
            match idx.strand_by_read.get(id) {
                Some(prokadiff_evidence::Strand::Plus) => ref_fwd += 1,
                _ => ref_rev += 1,
            }
        } else {
            other += 1;
        }
    }

    AlleleCounts::new(
        StrandCounts::new(ref_fwd, ref_rev),
        StrandCounts::new(alt_fwd, alt_rev),
        other,
    )
}

/// Per-read DEL allele classification at column `idx`, returning the 1 bp and
/// 2 bp allele counts in one pass.
///
/// M11 Round-4 (P1): deletion spans come from the structural `DeletionObs`
/// records the builder attaches at each run's START column after final
/// right-normalization. Spans are NEVER re-inferred from surviving `b'-'`
/// base observations: end-trimming legitimately removes edge dashes (e.g. a
/// `70M2D1M` read whose deletion right-normalizes onto the read's 3′ end
/// loses its last dash) and would otherwise fabricate a 1 bp allele from a
/// 2 bp event.
///
/// * ALT(1) = reads with a `DeletionObs` AT idx of span 1.
/// * ALT(2) = reads with a `DeletionObs` AT idx of span 2.
/// * OTHER  = base mismatches at idx; `b'-'` observations whose deletion
///   event is anchored upstream (continuations, counted as ALT at the run
///   start column instead); `b'-'` observations without a matching event
///   record at idx; deletion events with span >= 3 (outside the authorized
///   <=2 bp scope).
/// * Depth  = read-union formula:
///   `depth(x) = |{read_id : base obs at x} ∪ {read_id : deletion event covers x, no base obs at x}|`
///   A deletion event recorded at column `j` covers `j .. j + span - 1`
///   (recorded spans are 1–2, so only `j ∈ {idx-1, idx}` can cover idx).
///   Per-read_id union, no double counting; REF + ALT + OTHER == depth by
///   construction of the three passes below.
fn extract_del_counts(
    cols: &[prokadiff_evidence::PileupColumn],
    idx: usize,
) -> (AlleleCounts, AlleleCounts) {
    let col = &cols[idx];
    let ref_b = col.ref_base.to_ascii_uppercase();

    let mut one = AlleleCounts::default();
    let mut two = AlleleCounts::default();

    // Reads with a deletion event record anchored at idx (at most one per
    // read per column by builder construction).
    let event_at: std::collections::HashSet<u64> =
        col.deletions.iter().map(|ev| ev.read_id).collect();

    // Pass 1: every read with a base observation at idx, counted exactly once.
    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for obs in &col.observations {
        if !seen.insert(obs.read_id) {
            continue; // defensive: one base observation per read per column
        }
        let base = obs.base.to_ascii_uppercase();
        if base == b'-' {
            if event_at.contains(&obs.read_id) {
                // Counted from its event record in pass 2, never twice.
                continue;
            }
            // Continuation of an upstream event (or a dash observation
            // without a record): never an ALT at idx.
            one.other += 1;
            two.other += 1;
        } else if base == ref_b {
            let plus = obs.strand == prokadiff_evidence::Strand::Plus;
            if plus {
                one.ref_counts.fwd += 1;
                two.ref_counts.fwd += 1;
            } else {
                one.ref_counts.rev += 1;
                two.ref_counts.rev += 1;
            }
        } else {
            one.other += 1;
            two.other += 1;
        }
    }

    // Pass 2: deletion events anchored at idx (ALT / out-of-scope OTHER).
    // Strand comes from the record itself so a read whose base observation
    // was end-trimmed is still counted.
    for ev in &col.deletions {
        let plus = ev.strand == prokadiff_evidence::Strand::Plus;
        match ev.span {
            1 => {
                if plus {
                    one.alt_counts.fwd += 1;
                } else {
                    one.alt_counts.rev += 1;
                }
                // For the 2 bp candidate this read is neither REF nor ALT
                // (denominator symmetry).
                two.other += 1;
            }
            2 => {
                if plus {
                    two.alt_counts.fwd += 1;
                } else {
                    two.alt_counts.rev += 1;
                }
                one.other += 1;
            }
            // >= 3 bp event: outside the authorized subclonal scope.
            _ => {
                one.other += 1;
                two.other += 1;
            }
        }
    }

    // Pass 3: union-depth closure. Reads whose span>=2 event anchored at
    // idx-1 COVERS idx but that have no base observation here (their '-' at
    // idx was end-trimmed): depth OTHER, never ALT at idx.
    if idx > 0 {
        for ev in &cols[idx - 1].deletions {
            if ev.span < 2 || !seen.insert(ev.read_id) || event_at.contains(&ev.read_id) {
                continue;
            }
            one.other += 1;
            two.other += 1;
        }
    }

    (one, two)
}

/// Extracts candidate small variant observations directly from starter and edited pileups.
///
/// M11 Round-3 repair: allele eligibility is decided per read (via
/// [`ColumnReadIndex`]) before aggregation; deletion spans come from
/// per-read `DeletionObs` event records (M11 Round-4: adjacent disjoint
/// deletions stay separate 1 bp alleles, and end-trimmed edge dashes can
/// no longer shorten a 2 bp event); INS REF/ALT come from disjoint
/// per-read sets with no aggregate subtraction.
pub fn extract_subclonal_candidates_from_pileups(
    contig_names: &[String],
    starter_pileups: &[prokadiff_evidence::ContigPileup],
    edited_pileups: &[prokadiff_evidence::ContigPileup],
    min_af: f64,
    min_alt_reads: u64,
) -> Vec<RawSubclonalObservation> {
    let mut candidates = Vec::new();
    let n_contigs = contig_names.len().min(starter_pileups.len()).min(edited_pileups.len());

    for c_idx in 0..n_contigs {
        let contig = &contig_names[c_idx];
        let s_cols = &starter_pileups[c_idx].columns;
        let e_cols = &edited_pileups[c_idx].columns;
        let n_pos = s_cols.len().min(e_cols.len());

        for idx in 0..n_pos {
            let pos = idx as u64 + 1;
            let s_col = &s_cols[idx];
            let e_col = &e_cols[idx];
            let ref_char = (e_col.ref_base as char).to_ascii_uppercase();

            // Per-read indexes, built once per column and shared by the four
            // SNP candidates and every INS candidate at this position.
            let e_idx = ColumnReadIndex::new(e_col);
            let s_idx = ColumnReadIndex::new(s_col);

            // 1. Candidate SNPs
            for &alt_b in b"ACGT" {
                if alt_b == e_col.ref_base.to_ascii_uppercase() {
                    continue;
                }
                let e_counts =
                    extract_snp_counts_from_index(&e_idx, e_col.ref_base.to_ascii_uppercase(), alt_b);
                let e_alt = e_counts.alt_counts.total();
                let e_depth = e_counts.total_depth();
                let e_af = if e_depth > 0 { e_alt as f64 / e_depth as f64 } else { 0.0 };

                let s_counts =
                    extract_snp_counts_from_index(&s_idx, s_col.ref_base.to_ascii_uppercase(), alt_b);
                let s_alt = s_counts.alt_counts.total();
                let s_depth = s_counts.total_depth();
                let s_af = if s_depth > 0 { s_alt as f64 / s_depth as f64 } else { 0.0 };

                let is_e_cand = e_alt >= min_alt_reads && e_af >= min_af && e_af < 0.75;
                let is_s_cand = s_alt >= min_alt_reads && s_af >= min_af && s_af < 0.75;

                if is_e_cand || is_s_cand {
                    candidates.push(RawSubclonalObservation {
                        contig: contig.clone(),
                        pos,
                        ref_seq: ref_char.to_string(),
                        alt_seq: (alt_b as char).to_ascii_uppercase().to_string(),
                        variant_type: "SNP".to_string(),
                        starter_counts: s_counts,
                        edited_counts: e_counts,
                        starter_error_rate: calc_window_error(s_cols, idx),
                        edited_error_rate: calc_window_error(e_cols, idx),
                    });
                }
            }

            // 2. Candidate DEL (1 bp or 2 bp).
            //
            // Per-read deletion-event identity (P1-1, Round-4 P1): each
            // read's deletion is ONE event recorded by the builder at its
            // normalized start column; span is taken from the structural
            // `DeletionObs` record, never re-inferred from surviving `b'-'
            // base observations. Events starting at idx with span 1 and 2
            // are counted as separate alleles; span >= 3 is outside the
            // authorized subclonal scope and only feeds OTHER/depth.
            // Adjacent disjoint deletions (20 reads deleting only x plus 20
            // different reads deleting only x+1) stay two 1 bp candidates
            // instead of merging into a false 2 bp allele.
            let (e_del1, e_del2) = extract_del_counts(e_cols, idx);
            let (s_del1, s_del2) = extract_del_counts(s_cols, idx);

            for (del_len, e_counts, s_counts) in [
                (1usize, &e_del1, &s_del1),
                (2, &e_del2, &s_del2),
            ] {
                if del_len == 2 && idx + 1 >= n_pos {
                    continue;
                }
                let e_alt = e_counts.alt_counts.total();
                let e_depth = e_counts.total_depth();
                let e_af = if e_depth > 0 { e_alt as f64 / e_depth as f64 } else { 0.0 };

                let s_alt = s_counts.alt_counts.total();
                let s_depth = s_counts.total_depth();
                let s_af = if s_depth > 0 { s_alt as f64 / s_depth as f64 } else { 0.0 };

                let is_e_cand = e_alt >= min_alt_reads && e_af >= min_af && e_af < 0.75;
                let is_s_cand = s_alt >= min_alt_reads && s_af >= min_af && s_af < 0.75;

                if is_e_cand || is_s_cand {
                    let ref_seq = if del_len == 2 {
                        format!(
                            "{}{}",
                            ref_char,
                            (e_cols[idx + 1].ref_base as char).to_ascii_uppercase()
                        )
                    } else {
                        ref_char.to_string()
                    };
                    candidates.push(RawSubclonalObservation {
                        contig: contig.clone(),
                        pos,
                        ref_seq,
                        alt_seq: "-".to_string(),
                        variant_type: "DEL".to_string(),
                        starter_counts: *s_counts,
                        edited_counts: *e_counts,
                        starter_error_rate: calc_window_error(s_cols, idx),
                        edited_error_rate: calc_window_error(e_cols, idx),
                    });
                }
            }

            // 3. Candidate 1-2 bp INS
            let mut seen_ins = std::collections::BTreeSet::new();
            for ins in e_col.insertions.iter().chain(s_col.insertions.iter()) {
                if (1..=2).contains(&ins.seq.len()) {
                    seen_ins.insert(ins.seq.to_ascii_uppercase());
                }
            }

            for ins_seq in seen_ins {
                let e_counts = extract_ins_counts(&e_idx, e_col.ref_base.to_ascii_uppercase(), &ins_seq);
                let e_alt = e_counts.alt_counts.total();
                let e_depth = e_counts.total_depth();
                let e_af = if e_depth > 0 { e_alt as f64 / e_depth as f64 } else { 0.0 };

                let s_counts = extract_ins_counts(&s_idx, s_col.ref_base.to_ascii_uppercase(), &ins_seq);
                let s_alt = s_counts.alt_counts.total();
                let s_depth = s_counts.total_depth();
                let s_af = if s_depth > 0 { s_alt as f64 / s_depth as f64 } else { 0.0 };

                let is_e_cand = e_alt >= min_alt_reads && e_af >= min_af && e_af < 0.75;
                let is_s_cand = s_alt >= min_alt_reads && s_af >= min_af && s_af < 0.75;

                if is_e_cand || is_s_cand {
                    candidates.push(RawSubclonalObservation {
                        contig: contig.clone(),
                        pos,
                        ref_seq: "-".to_string(),
                        alt_seq: String::from_utf8_lossy(&ins_seq).into_owned(),
                        variant_type: "INS".to_string(),
                        starter_counts: s_counts,
                        edited_counts: e_counts,
                        starter_error_rate: calc_window_error(s_cols, idx),
                        edited_error_rate: calc_window_error(e_cols, idx),
                    });
                }
            }
        }
    }

    candidates
}

#[cfg(test)]
mod tests {
    use super::*;
    use prokadiff_evidence::subclonal::StrandCounts;

    #[test]
    fn test_evaluate_subclonal_candidates_empty() {
        let params = SubclonalParams::default();
        let res = evaluate_subclonal_candidates(Vec::new(), &params);
        assert!(res.is_empty());
    }

    #[test]
    fn test_evaluate_subclonal_candidates_single_edited_specific() {
        let params = SubclonalParams::default();
        let obs = RawSubclonalObservation {
            contig: "chr1".to_string(),
            pos: 100,
            ref_seq: "A".to_string(),
            alt_seq: "G".to_string(),
            variant_type: "SNP".to_string(),
            starter_counts: AlleleCounts::new(
                StrandCounts::new(50, 50),
                StrandCounts::new(0, 0),
                0,
            ),
            edited_counts: AlleleCounts::new(
                StrandCounts::new(45, 45),
                StrandCounts::new(5, 5),
                0,
            ),
            starter_error_rate: 0.005,
            edited_error_rate: 0.005,
        };

        let res = evaluate_subclonal_candidates(vec![obs], &params);
        assert_eq!(res.len(), 1);
        let v = &res[0];
        assert_eq!(v.classification, SubclonalClass::EditedSpecific);
        assert_eq!(v.starter_af, Some(0.0));
        assert_eq!(v.edited_af, Some(0.10));
        assert!(v.delta_ci_low > 0.0);
        assert!(v.starter_exclusion_bound.is_some());
        assert!(!v.flags.contains(&"STRAND_BIAS".to_string()));
    }

    // ---- M11 Round-3 repair: per-read allele classification (P1-1 / P1-2) ----

    use prokadiff_evidence::Strand as RaStrand;
    use prokadiff_evidence::{BaseObs, ContigPileup, DeletionObs, InsertionObs, PileupColumn};

    fn ref_cols(ref_seq: &[u8]) -> Vec<PileupColumn> {
        ref_seq
            .iter()
            .map(|&b| PileupColumn {
                ref_base: b,
                observations: Vec::new(),
                insertions: Vec::new(),
                deletions: Vec::new(),
            })
            .collect()
    }

    fn add_base(
        cols: &mut [PileupColumn],
        idx: usize,
        read_id: u64,
        base: u8,
        strand: RaStrand,
    ) {
        cols[idx].observations.push(BaseObs {
            base,
            strand,
            read_id,
        });
    }

    fn add_ins(cols: &mut [PileupColumn], idx: usize, read_id: u64, seq: &[u8], strand: RaStrand) {
        cols[idx].insertions.push(InsertionObs {
            seq: seq.to_vec(),
            strand,
            read_id,
        });
    }

    /// One read's L-bp deletion event starting at `idx`, exactly as the real
    /// builder emits it: dash `BaseObs` at every covered column PLUS one
    /// `DeletionObs` record (span = L) at the run's start column.
    fn add_del(
        cols: &mut [PileupColumn],
        idx: usize,
        read_id: u64,
        span: usize,
        strand: RaStrand,
    ) {
        for d in 0..span {
            cols[idx + d].observations.push(BaseObs {
                base: b'-',
                strand,
                read_id,
            });
        }
        cols[idx].deletions.push(DeletionObs {
            span: span as u32,
            strand,
            read_id,
        });
    }

    fn pileup(cols: Vec<PileupColumn>) -> ContigPileup {
        let n = cols.len();
        ContigPileup {
            columns: cols,
            unique_depth: vec![0; n],
            total_depth: vec![0; n],
            splits: Vec::new(),
        }
    }

    fn extract(
        s_cols: Vec<PileupColumn>,
        e_cols: Vec<PileupColumn>,
    ) -> Vec<RawSubclonalObservation> {
        extract_subclonal_candidates_from_pileups(
            &["chr1".to_string()],
            &[pileup(s_cols)],
            &[pileup(e_cols)],
            0.05,
            3,
        )
    }

    fn find_candidate<'a>(
        cands: &'a [RawSubclonalObservation],
        pos: u64,
        ref_seq: &str,
        alt_seq: &str,
    ) -> &'a RawSubclonalObservation {
        cands
            .iter()
            .find(|c| c.pos == pos && c.ref_seq == ref_seq && c.alt_seq == alt_seq)
            .unwrap_or_else(|| panic!("missing candidate at pos {pos} ref {ref_seq} alt {alt_seq}; got {cands:?}"))
    }

    const STRETCH: &[u8] = b"ACGTACGTACGT"; // index 4 = 'A', index 5 = 'C'

    /// P1-1: reads deleting only x plus different reads deleting only x+1 are
    /// two separate 1 bp alleles, never a merged false 2 bp deletion.
    #[test]
    fn disjoint_adjacent_deletions_stay_separate_one_bp_alleles() {
        let s_cols = ref_cols(STRETCH);
        let mut e_cols = ref_cols(STRETCH);

        // Reads 0..10 delete only reference index 4 ('A'): '-' at 4, 'C' at 5.
        for id in 0..10u64 {
            let strand = if id % 2 == 0 { RaStrand::Plus } else { RaStrand::Minus };
            add_del(&mut e_cols, 4, id, 1, strand);
            add_base(&mut e_cols, 5, id, b'C', strand);
        }
        // Reads 10..20 delete only reference index 5 ('C'): 'A' at 4, '-' at 5.
        for id in 10..20u64 {
            let strand = if id % 2 == 0 { RaStrand::Plus } else { RaStrand::Minus };
            add_base(&mut e_cols, 4, id, b'A', strand);
            add_del(&mut e_cols, 5, id, 1, strand);
        }
        // 20 pure-REF reads.
        for id in 20..40u64 {
            add_base(&mut e_cols, 4, id, b'A', RaStrand::Plus);
            add_base(&mut e_cols, 5, id, b'C', RaStrand::Plus);
        }
        // Starter: 40 pure-REF reads (ids 100..140 never collide with edited).
        let mut s_cols = s_cols;
        for id in 100..140u64 {
            add_base(&mut s_cols, 4, id, b'A', RaStrand::Plus);
            add_base(&mut s_cols, 5, id, b'C', RaStrand::Plus);
        }

        let cands = extract(s_cols, e_cols);

        // Exactly two 1 bp DEL candidates; no 2 bp candidate anywhere.
        let del1_a = find_candidate(&cands, 5, "A", "-");
        let del1_c = find_candidate(&cands, 6, "C", "-");
        assert_eq!(del1_a.variant_type, "DEL");
        assert_eq!(del1_c.variant_type, "DEL");
        assert_eq!(del1_a.edited_counts.alt_counts.total(), 10);
        assert_eq!(del1_c.edited_counts.alt_counts.total(), 10);
        // Denominator: all reads with any observation at the column (40).
        assert_eq!(del1_a.edited_counts.total_depth(), 40);
        assert_eq!(del1_c.edited_counts.total_depth(), 40);
        assert_eq!(del1_a.starter_counts.total_depth(), 40);
        assert_eq!(del1_a.starter_counts.alt_counts.total(), 0);
        assert!(
            !cands.iter().any(|c| c.variant_type == "DEL" && c.ref_seq.len() == 2),
            "adjacent disjoint deletions must not produce a 2 bp candidate: {cands:?}"
        );
    }

    /// P1-1: span-1 and span-2 deletion runs starting at the same column keep
    /// separate allele identities and correct denominators.
    #[test]
    fn mixed_deletion_lengths_at_one_locus_separate_alleles() {
        let s_cols = ref_cols(STRETCH);
        let mut e_cols = ref_cols(STRETCH);

        // 8 reads delete exactly 1 bp at index 4.
        for id in 0..8u64 {
            let strand = if id % 2 == 0 { RaStrand::Plus } else { RaStrand::Minus };
            add_del(&mut e_cols, 4, id, 1, strand);
            add_base(&mut e_cols, 5, id, b'C', strand);
        }
        // 12 reads delete exactly 2 bp (indices 4 and 5).
        for id in 8..20u64 {
            let strand = if id % 2 == 0 { RaStrand::Plus } else { RaStrand::Minus };
            add_del(&mut e_cols, 4, id, 2, strand);
        }
        // 20 pure-REF reads.
        for id in 20..40u64 {
            add_base(&mut e_cols, 4, id, b'A', RaStrand::Plus);
            add_base(&mut e_cols, 5, id, b'C', RaStrand::Plus);
        }
        let mut s_cols = s_cols;
        for id in 100..140u64 {
            add_base(&mut s_cols, 4, id, b'A', RaStrand::Plus);
            add_base(&mut s_cols, 5, id, b'C', RaStrand::Plus);
        }

        let cands = extract(s_cols, e_cols);

        let del1 = find_candidate(&cands, 5, "A", "-");
        let del2 = find_candidate(&cands, 5, "AC", "-");
        assert_eq!(del1.edited_counts.alt_counts.total(), 8);
        assert_eq!(del2.edited_counts.alt_counts.total(), 12);
        // Shared column depth denominator for both alleles.
        assert_eq!(del1.edited_counts.total_depth(), 40);
        assert_eq!(del2.edited_counts.total_depth(), 40);
        assert_eq!(del1.starter_counts.alt_counts.total(), 0);
        assert_eq!(del2.starter_counts.alt_counts.total(), 0);
        // The 2 bp reads must NOT surface as a false 1 bp allele at pos 6:
        // their run starts at index 4, so index 5 is a continuation.
        assert!(
            !cands.iter().any(|c| c.pos == 6 && c.variant_type == "DEL"),
            "continuation of a 2 bp run must not become a 1 bp candidate at pos 6: {cands:?}"
        );
    }

    /// P1-2: an insertion whose anchor base was repeat-masked (no base
    /// observation) is still ALT via its insertion entry, and REF is a strict
    /// set difference — no saturating_sub stealing from pure-REF reads.
    #[test]
    fn insertion_with_masked_anchor_counts_alt_without_ref_stealing() {
        let s_cols = ref_cols(STRETCH);
        let mut e_cols = ref_cols(STRETCH);

        // 10 insertion-only reads: entry "GT" at index 4, NO base observation.
        for id in 0..10u64 {
            let strand = if id % 2 == 0 { RaStrand::Plus } else { RaStrand::Minus };
            add_ins(&mut e_cols, 4, id, b"GT", strand);
        }
        // 30 pure-REF reads with base observations.
        for id in 10..40u64 {
            add_base(&mut e_cols, 4, id, b'A', RaStrand::Plus);
        }
        let mut s_cols = s_cols;
        for id in 100..140u64 {
            add_base(&mut s_cols, 4, id, b'A', RaStrand::Plus);
        }

        let cands = extract(s_cols, e_cols);
        let ins = find_candidate(&cands, 5, "-", "GT");
        assert_eq!(ins.variant_type, "INS");
        assert_eq!(ins.edited_counts.alt_counts.total(), 10);
        // REF must be exactly the 30 clean reads — none stolen by insertions.
        assert_eq!(ins.edited_counts.ref_counts.total(), 30);
        // Depth counts insertion-only reads once via their entry: 30 + 10 = 40.
        assert_eq!(ins.edited_counts.total_depth(), 40);
        assert_eq!(ins.starter_counts.total_depth(), 40);
        assert_eq!(ins.starter_counts.alt_counts.total(), 0);
        // The union invariant: depth == |base-obs reads ∪ insertion reads|.
        let e_cols_union = {
            let mut cols = ref_cols(STRETCH);
            for id in 0..10u64 {
                add_ins(&mut cols, 4, id, b"GT", RaStrand::Plus);
            }
            for id in 10..40u64 {
                add_base(&mut cols, 4, id, b'A', RaStrand::Plus);
            }
            cols
        };
        assert_eq!(
            ColumnReadIndex::new(&e_cols_union[4]).depth_union(),
            40,
            "masked-anchor reads count once via their insertion entry"
        );
    }

    /// P1-2: two different insertion oligos at the same anchor each count as
    /// ALT for their own candidate and OTHER for the rival; no REF stealing.
    #[test]
    fn two_insertion_sequences_count_separately_without_ref_stealing() {
        let s_cols = ref_cols(STRETCH);
        let mut e_cols = ref_cols(STRETCH);

        for id in 0..8u64 {
            add_ins(&mut e_cols, 4, id, b"GT", RaStrand::Plus);
        }
        for id in 8..16u64 {
            add_ins(&mut e_cols, 4, id, b"AC", RaStrand::Minus);
        }
        for id in 16..40u64 {
            add_base(&mut e_cols, 4, id, b'A', RaStrand::Plus);
        }
        let mut s_cols = s_cols;
        for id in 100..140u64 {
            add_base(&mut s_cols, 4, id, b'A', RaStrand::Plus);
        }

        let cands = extract(s_cols, e_cols);
        let ins_gt = find_candidate(&cands, 5, "-", "GT");
        let ins_ac = find_candidate(&cands, 5, "-", "AC");

        assert_eq!(ins_gt.edited_counts.alt_counts.total(), 8);
        assert_eq!(ins_ac.edited_counts.alt_counts.total(), 8);
        // Each candidate sees the rival oligo reads as OTHER, not REF.
        assert_eq!(ins_gt.edited_counts.other, 8);
        assert_eq!(ins_ac.edited_counts.other, 8);
        assert_eq!(ins_gt.edited_counts.ref_counts.total(), 24);
        assert_eq!(ins_ac.edited_counts.ref_counts.total(), 24);
        assert_eq!(ins_gt.edited_counts.total_depth(), 40);
        assert_eq!(ins_ac.edited_counts.total_depth(), 40);
    }

    /// Deletion runs of length >= 3 are outside the authorized <=2 bp scope:
    /// they feed OTHER/depth only and never produce <=2 bp candidates. The
    /// real builder diverts >2 bp CIGAR deletions to the JC path (no dashes,
    /// no event records); this test keeps the defensive fallback for dash
    /// observations without a matching `DeletionObs` record.
    #[test]
    fn three_bp_deletion_reads_stay_in_other_without_candidates() {
        let s_cols = ref_cols(STRETCH);
        let mut e_cols = ref_cols(STRETCH);

        // 10 reads with a 3 bp deletion at indices 4..7.
        for id in 0..10u64 {
            let strand = if id % 2 == 0 { RaStrand::Plus } else { RaStrand::Minus };
            add_base(&mut e_cols, 4, id, b'-', strand);
            add_base(&mut e_cols, 5, id, b'-', strand);
            add_base(&mut e_cols, 6, id, b'-', strand);
        }
        for id in 10..40u64 {
            add_base(&mut e_cols, 4, id, b'A', RaStrand::Plus);
            add_base(&mut e_cols, 5, id, b'C', RaStrand::Plus);
            add_base(&mut e_cols, 6, id, b'G', RaStrand::Plus);
        }
        let mut s_cols = s_cols;
        for id in 100..140u64 {
            add_base(&mut s_cols, 4, id, b'A', RaStrand::Plus);
            add_base(&mut s_cols, 5, id, b'C', RaStrand::Plus);
            add_base(&mut s_cols, 6, id, b'G', RaStrand::Plus);
        }

        let cands = extract(s_cols, e_cols.clone());
        assert!(
            !cands.iter().any(|c| c.variant_type == "DEL"),
            "3 bp deletion reads must not produce <=2 bp candidates: {cands:?}"
        );

        // The out-of-scope reads are counted into OTHER with correct depth.
        let (del1, del2) = extract_del_counts(&e_cols, 4);
        assert_eq!(del1.alt_counts.total(), 0);
        assert_eq!(del2.alt_counts.total(), 0);
        assert_eq!(del1.other, 10);
        assert_eq!(del1.ref_counts.total(), 30);
        assert_eq!(del1.total_depth(), 40);
        assert_eq!(del2.total_depth(), 40);
    }

    /// SNP eligibility is insertion-aware: a read carrying an insertion entry
    /// at the column is OTHER for the SNP, never ALT or REF.
    #[test]
    fn snp_excludes_insertion_carrier_reads() {
        let mut col = ref_cols(STRETCH).remove(4);
        // 10 reads: base == ALT 'G' but with an insertion entry at this column.
        for id in 0..10u64 {
            col.observations.push(BaseObs {
                base: b'G',
                strand: RaStrand::Plus,
                read_id: id,
            });
            col.insertions.push(InsertionObs {
                seq: b"T".to_vec(),
                strand: RaStrand::Plus,
                read_id: id,
            });
        }
        for id in 10..40u64 {
            col.observations.push(BaseObs {
                base: b'A',
                strand: RaStrand::Plus,
                read_id: id,
            });
        }

        let idx = ColumnReadIndex::new(&col);
        let counts = extract_snp_counts_from_index(&idx, col.ref_base.to_ascii_uppercase(), b'G');
        assert_eq!(counts.alt_counts.total(), 0, "insertion carriers are not SNP ALT");
        assert_eq!(counts.ref_counts.total(), 30);
        assert_eq!(counts.other, 10);
        assert_eq!(counts.total_depth(), 40);
    }

    /// End-to-end through the real pileup builder: per-read identity assigned
    /// by the builder lets disjoint adjacent deletions stay separate alleles.
    #[test]
    fn builder_assigned_read_ids_keep_adjacent_deletions_disjoint() {
        use prokadiff_evidence::pileup::{apply_read, AlignedRead, CigarKind, CigarOp};

        let ref_seq = STRETCH.to_vec();
        let mut e_cols = ref_cols(&ref_seq);
        let mut unique_depth = vec![0u32; ref_seq.len()];
        let mut total_depth = vec![0u32; ref_seq.len()];
        let mut splits = Vec::new();
        let mut clips = Vec::new();

        // 10 reads deleting only reference index 4 (1D at position 5).
        for ordinal in 0..10u64 {
            let read = AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus: ordinal % 2 == 1,
                seq: b"ACGTCGTACGT".to_vec(),
                cigar: vec![
                    CigarOp { kind: CigarKind::Match, len: 4 },
                    CigarOp { kind: CigarKind::Del, len: 1 },
                    CigarOp { kind: CigarKind::Match, len: 7 },
                ],
                mapq: 40,
                molecule_id: 0, // deliberately identical: builder ordinal is the identity
            };
            apply_read(
                &read,
                ordinal,
                &mut e_cols,
                &mut unique_depth,
                &mut total_depth,
                &mut splits,
                &mut clips,
                Some(&ref_seq),
            );
        }
        // 10 reads deleting only reference index 5 (1D at position 6).
        for ordinal in 10..20u64 {
            let read = AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus: ordinal % 2 == 1,
                seq: b"ACGTACTACGT".to_vec(),
                cigar: vec![
                    CigarOp { kind: CigarKind::Match, len: 5 },
                    CigarOp { kind: CigarKind::Del, len: 1 },
                    CigarOp { kind: CigarKind::Match, len: 6 },
                ],
                mapq: 40,
                molecule_id: 0,
            };
            apply_read(
                &read,
                ordinal,
                &mut e_cols,
                &mut unique_depth,
                &mut total_depth,
                &mut splits,
                &mut clips,
                Some(&ref_seq),
            );
        }

        let mut s_cols = ref_cols(&ref_seq);
        for ordinal in 100..140u64 {
            let read = AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus: false,
                seq: ref_seq.clone(),
                cigar: vec![CigarOp { kind: CigarKind::Match, len: 12 }],
                mapq: 40,
                molecule_id: 0,
            };
            apply_read(
                &read,
                ordinal,
                &mut s_cols,
                &mut vec![0u32; ref_seq.len()],
                &mut vec![0u32; ref_seq.len()],
                &mut Vec::new(),
                &mut Vec::new(),
                Some(&ref_seq),
            );
        }

        let cands = extract(s_cols, e_cols);
        assert!(
            !cands.iter().any(|c| c.variant_type == "DEL" && c.ref_seq.len() == 2),
            "builder path must not merge disjoint adjacent deletions: {cands:?}"
        );
        let del_a = find_candidate(&cands, 5, "A", "-");
        let del_c = find_candidate(&cands, 6, "C", "-");
        assert_eq!(del_a.edited_counts.alt_counts.total(), 10);
        assert_eq!(del_c.edited_counts.alt_counts.total(), 10);
    }

    /// Debug invariant (M11 Round-4 P1) that replaces the removed dash-run
    /// inference: every surviving `b'-'` base observation must belong to a
    /// `DeletionObs` event recorded at some column j <= idx with
    /// j + span > idx (the event covering it). Builder trimming may remove
    /// dashes, never events.
    #[cfg(test)]
    fn assert_every_dash_is_covered_by_an_event(cols: &[PileupColumn]) {
        for (idx, col) in cols.iter().enumerate() {
            for obs in &col.observations {
                if obs.base != b'-' {
                    continue;
                }
                let covered = (0..=idx).any(|j| {
                    cols[j].deletions.iter().any(|ev| {
                        ev.read_id == obs.read_id && j + ev.span as usize > idx
                    })
                });
                assert!(
                    covered,
                    "dash of read {} at column {idx} has no covering DeletionObs",
                    obs.read_id
                );
            }
        }
    }

    /// M11 Round-4 (P1) REQUIRED regression: the Codex-Sol R4 counterexample
    /// through the real builder.
    ///
    /// Reference `P + "ACA"` (P = 70 non-homopolymer bases, 73 bp total), 10
    /// edited reads `P + "A"` with CIGAR `70M2D1M`. The 2 bp deletion
    /// right-normalizes onto refs 71..72, the trailing dash lands on the
    /// read's last aligned base and is end-trimmed, so only ONE `'-'`
    /// observation survives per deletion read. The old dash-run inference
    /// saw that 1-dash run and emitted a false 1 bp DEL; the `DeletionObs`
    /// event record preserves span 2 at the normalized start, so the
    /// classify output is exactly ONE 2 bp DEL candidate.
    #[test]
    fn trimmed_edge_deletion_stays_a_two_bp_candidate() {
        use prokadiff_evidence::pileup::{apply_read, AlignedRead, CigarKind, CigarOp};

        let mut ref_seq = Vec::new();
        for i in 0..70 {
            ref_seq.push(b"ACGT"[i % 4]);
        }
        ref_seq.extend_from_slice(b"ACA"); // refs 70='A', 71='C', 72='A'
        let mut e_cols = ref_cols(&ref_seq);
        let mut unique_depth = vec![0u32; ref_seq.len()];
        let mut total_depth = vec![0u32; ref_seq.len()];
        let mut splits = Vec::new();
        let mut clips = Vec::new();

        let mut read_seq = ref_seq[..70].to_vec();
        read_seq.push(b'A'); // 71 bp: P + "A"
        for ordinal in 0..10u64 {
            let read = AlignedRead {
                contig_idx: 0,
                ref_start_0: 0,
                minus: ordinal % 2 == 1,
                seq: read_seq.clone(),
                cigar: vec![
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 70,
                    },
                    CigarOp {
                        kind: CigarKind::Del,
                        len: 2,
                    },
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 1,
                    },
                ],
                mapq: 40,
                molecule_id: 0,
            };
            apply_read(
                &read,
                ordinal,
                &mut e_cols,
                &mut unique_depth,
                &mut total_depth,
                &mut splits,
                &mut clips,
                Some(&ref_seq),
            );
        }
        // 30 pure-REF edited reads + 40 pure-REF starter reads.
        let mut s_cols = ref_cols(&ref_seq);
        for (cols, range) in [
            (&mut e_cols, 10..40u64),
            (&mut s_cols, 100..140u64),
        ] {
            for ordinal in range {
                let read = AlignedRead {
                    contig_idx: 0,
                    ref_start_0: 0,
                    minus: false,
                    seq: ref_seq.clone(),
                    cigar: vec![CigarOp {
                        kind: CigarKind::Match,
                        len: ref_seq.len(),
                    }],
                    mapq: 40,
                    molecule_id: 0,
                };
                apply_read(
                    &read,
                    ordinal,
                    cols,
                    &mut vec![0u32; ref_seq.len()],
                    &mut vec![0u32; ref_seq.len()],
                    &mut Vec::new(),
                    &mut Vec::new(),
                    Some(&ref_seq),
                );
            }
        }

        // Builder-level assertions: the trim really removed each deletion's
        // final dash, and each read kept one span-2 event at the normalized
        // start ref 71.
        for id in 0..10u64 {
            let dashes: Vec<usize> = e_cols
                .iter()
                .enumerate()
                .filter(|(_, c)| c.observations.iter().any(|o| o.read_id == id && o.base == b'-'))
                .map(|(i, _)| i)
                .collect();
            assert_eq!(dashes, vec![71], "read {id}: exactly one surviving dash");
            assert_eq!(e_cols[71].deletions.iter().filter(|ev| ev.read_id == id).count(), 1);
            assert_eq!(e_cols[71].deletions.iter().find(|ev| ev.read_id == id).unwrap().span, 2);
        }
        assert_eq!(e_cols[71].deletions.len(), 10, "one event per read per run");
        assert_every_dash_is_covered_by_an_event(&e_cols);

        // Classify: exactly one 2 bp DEL candidate, no false 1 bp allele.
        let cands = extract(s_cols, e_cols);
        let dels: Vec<&RawSubclonalObservation> =
            cands.iter().filter(|c| c.variant_type == "DEL").collect();
        assert_eq!(dels.len(), 1, "exactly one DEL candidate: {cands:?}");
        let del = dels[0];
        assert_eq!(del.pos, 72, "normalized start, 1-based");
        assert_eq!(del.ref_seq, "CA");
        assert_eq!(del.alt_seq, "-");
        assert_eq!(del.edited_counts.alt_counts.total(), 10);
        assert_eq!(del.edited_counts.total_depth(), 40);
        assert_eq!(del.starter_counts.alt_counts.total(), 0);

        // Through the statistical evaluation: exactly one published variant.
        let params = SubclonalParams::default();
        let res = evaluate_subclonal_candidates(cands, &params);
        assert_eq!(res.len(), 1, "exactly one published subclonal variant");
        assert_eq!(res[0].variant_type, "DEL");
        assert_eq!(res[0].ref_seq, "CA");
        assert_eq!(res[0].pos, 72);
        assert_eq!(res[0].classification, SubclonalClass::EditedSpecific);
    }

    /// M11 Round-4 (P1): union-depth formula. A read whose span-2 event was
    /// anchored at idx-1 and whose '-' at idx was end-trimmed away (deletion
    /// at the read's 3′ edge) still enters the depth at idx as OTHER — via
    /// its event coverage, never double counted.
    #[test]
    fn trimmed_event_coverage_enters_depth_as_other() {
        let s_cols = ref_cols(STRETCH);
        let mut e_cols = ref_cols(STRETCH);

        // Synthetic builder-faithful shape: 10 reads with a 2 bp event at
        // index 4 whose dash at index 5 was end-trimmed away (only the dash
        // at 4 plus the base observation at 3 survive — mirrors a deletion
        // run ending exactly at the read's last aligned base).
        for id in 0..10u64 {
            let strand = if id % 2 == 0 { RaStrand::Plus } else { RaStrand::Minus };
            add_base(&mut e_cols, 3, id, b'T', strand);
            add_base(&mut e_cols, 4, id, b'-', strand);
            e_cols[4].deletions.push(DeletionObs {
                span: 2,
                strand,
                read_id: id,
            });
        }
        // 30 pure-REF reads.
        for id in 10..40u64 {
            add_base(&mut e_cols, 4, id, b'A', RaStrand::Plus);
            add_base(&mut e_cols, 5, id, b'C', RaStrand::Plus);
        }
        let mut s_cols = s_cols;
        for id in 100..140u64 {
            add_base(&mut s_cols, 4, id, b'A', RaStrand::Plus);
            add_base(&mut s_cols, 5, id, b'C', RaStrand::Plus);
        }

        // At idx 5: the 10 reads have NO base observation but their span-2
        // event from idx 4 covers 5 -> union depth = 30 + 10 = 40, all OTHER.
        let (del1, del2) = extract_del_counts(&e_cols, 5);
        assert_eq!(del1.alt_counts.total(), 0);
        assert_eq!(del2.alt_counts.total(), 0);
        assert_eq!(del1.other, 10);
        assert_eq!(del2.other, 10);
        assert_eq!(del1.ref_counts.total(), 30);
        assert_eq!(del2.ref_counts.total(), 30);
        assert_eq!(del1.total_depth(), 40, "event-covered reads enter the union depth");
        assert_eq!(del2.total_depth(), 40);
        assert_eq!(del1.total_depth(), del2.total_depth());
    }
}
