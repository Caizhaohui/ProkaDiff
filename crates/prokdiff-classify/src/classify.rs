use prokdiff_gd::{GdEntry, GdKind, GenomeDiff};

use crate::homolog::{scan_homologs, HomologSite};
use crate::intended::{entry_intervals, mask_intended, IntendedEdit};
use crate::{is_product_mutation, is_structural, EditorKind, RefContig};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MutationClass {
    Structural,
    NearHomolog,
    ScatteredSnv,
}

impl MutationClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Structural => "structural",
            Self::NearHomolog => "near_homolog",
            Self::ScatteredSnv => "scattered_snv",
        }
    }
}

#[derive(Clone, Debug)]
pub struct ClassifyOptions {
    pub editor: EditorKind,
    pub spacer: Option<String>,
    pub pam: Option<String>,
    pub near_distance: u64,
    pub max_mismatches: u32,
    pub hypothesis: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClassifiedMutation {
    pub entry: GdEntry,
    pub class: MutationClass,
    pub pam_profile: Option<String>,
    pub offtarget_mismatch: Option<u32>,
    pub distance_to_site: Option<u64>,
    pub hypothesis: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct ClassifyResult {
    pub unintended: Vec<ClassifiedMutation>,
    pub intended_observed: Vec<GdEntry>,
    /// Number of rows in the `--intended` table (`intended.len()`). Zero when omitted.
    pub intended_declared: usize,
    pub starter_vs_ref: usize,
}

pub fn classify(
    edited: &GenomeDiff,
    starter: &GenomeDiff,
    intended: &[IntendedEdit],
    refs: &[RefContig],
    opts: &ClassifyOptions,
) -> ClassifyResult {
    let starter_muts: Vec<GdEntry> = starter
        .entries
        .iter()
        .filter(|e| is_product_mutation(e.kind))
        .cloned()
        .collect();
    let edited_only = GenomeDiff {
        metadata: edited.metadata.clone(),
        entries: edited
            .entries
            .iter()
            .filter(|e| is_product_mutation(e.kind))
            .cloned()
            .collect(),
    };
    let starter_only = GenomeDiff {
        metadata: starter.metadata.clone(),
        entries: starter_muts.clone(),
    };
    let diff = edited_only.subtract(&starter_only);
    let (remain, observed) = mask_intended(&diff.entries, intended);

    let pam_used = resolved_pam(opts);
    let sites = match (opts.editor, opts.spacer.as_deref(), pam_used.as_deref()) {
        (EditorKind::Dsb, _, _) | (_, None | Some(""), _) | (_, _, None | Some("")) => Vec::new(),
        (editor, Some(spacer), Some(pam)) => {
            scan_homologs(refs, spacer, pam, editor, opts.max_mismatches)
        }
    };

    let mut unintended = Vec::with_capacity(remain.len());
    for e in remain {
        unintended.push(label_one(e, &sites, opts, pam_used.as_deref()));
    }

    ClassifyResult {
        unintended,
        intended_observed: observed.into_iter().cloned().collect(),
        intended_declared: intended.len(),
        starter_vs_ref: starter_muts.len(),
    }
}

fn resolved_pam(opts: &ClassifyOptions) -> Option<String> {
    if let Some(p) = opts.pam.as_deref() {
        if !p.is_empty() {
            return Some(p.to_ascii_uppercase());
        }
    }
    opts.editor.default_pam().map(str::to_string)
}

fn label_one(
    e: &GdEntry,
    sites: &[HomologSite],
    opts: &ClassifyOptions,
    pam_used: Option<&str>,
) -> ClassifiedMutation {
    let del_size = if e.kind == GdKind::Del {
        e.fields.get(2).and_then(|s| s.parse().ok())
    } else {
        None
    };
    let (class, pam_profile, offtarget_mismatch, distance_to_site) =
        if is_structural(e.kind, del_size) {
            (MutationClass::Structural, None, None, None)
        } else if opts.editor != EditorKind::Dsb {
            if let Some((dist, mm)) = nearest_site(e, sites, opts.near_distance) {
                (
                    MutationClass::NearHomolog,
                    pam_used.map(str::to_string),
                    Some(mm),
                    Some(dist),
                )
            } else {
                (MutationClass::ScatteredSnv, None, None, None)
            }
        } else {
            (MutationClass::ScatteredSnv, None, None, None)
        };
    let hypothesis = if opts.hypothesis && class == MutationClass::ScatteredSnv {
        Some("sos_widney2014".into())
    } else {
        None
    };
    ClassifiedMutation {
        entry: e.clone(),
        class,
        pam_profile,
        offtarget_mismatch,
        distance_to_site,
        hypothesis,
    }
}

fn nearest_site(e: &GdEntry, sites: &[HomologSite], max_d: u64) -> Option<(u64, u32)> {
    let mut best: Option<(u64, u32)> = None;
    for (sid, a, b) in entry_intervals(e) {
        for s in sites {
            if s.seq_id != sid {
                continue;
            }
            let d = interval_distance(a, b, s.start, s.end);
            if d > max_d {
                continue;
            }
            match best {
                None => best = Some((d, s.mismatches)),
                Some((bd, bm)) if d < bd || (d == bd && s.mismatches < bm) => {
                    best = Some((d, s.mismatches));
                }
                _ => {}
            }
        }
    }
    best
}

fn interval_distance(a0: u64, a1: u64, b0: u64, b1: u64) -> u64 {
    if a0 <= b1 && b0 <= a1 {
        0
    } else if a0 > b1 {
        a0 - b1
    } else {
        b0 - a1
    }
}
