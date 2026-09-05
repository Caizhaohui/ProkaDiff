//! External Bowtie2 wrapper. Converts SAM → coordinate-sorted BAM with noodles.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as _;

use crate::error::{EvidenceError, Result};
use crate::jc_seq::MIN_CLIP_FOR_SEED;

/// Short seed for unmatched reads vs the reference (primary Stage 2).
/// Bowtie2 default `-L` is 20; unmatched 36 bp junction-supporting reads need
/// a hypersensitive pass. Observed oracle unmatched pass uses a single-digit
/// seed; we pin 9 here (do not copy breseq source).
pub const UNMATCHED_SEED_LEN: usize = 9;

/// Seed for candidate-junction second pass. Junction flanks are ~read_len/2
/// (~18 bp for Clonal 36 bp reads); default `-L 20` cannot seed those flanks
/// and spanning `min_side` collapses to ~1–3 bp (`min_side >= 14` → 0).
pub const JUNCTION_SEED_LEN: usize = 10;

/// Which Bowtie2 pipeline to run after `bowtie2-build`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignKind {
    /// Genome: Stage 1 `--local`, then unmatched re-aligned with [`UNMATCHED_SEED_LEN`].
    Primary,
    /// Junction constructs: [`JUNCTION_SEED_LEN`] so short flanks can seed.
    Junction,
}

/// One Bowtie2 argv (unit-tested; no subprocess).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AlignPass {
    PrimaryStage1,
    PrimaryStage2,
    Junction,
}

#[derive(Clone, Debug)]
pub struct FastqInput {
    pub files: Vec<PathBuf>,
}

impl FastqInput {
    pub fn se(path: impl Into<PathBuf>) -> Self {
        Self {
            files: vec![path.into()],
        }
    }

    pub fn pe(r1: impl Into<PathBuf>, r2: impl Into<PathBuf>) -> Self {
        Self {
            files: vec![r1.into(), r2.into()],
        }
    }
}

pub fn which_or_err(name: &'static str) -> Result<PathBuf> {
    which(name).ok_or(EvidenceError::MissingBinary(name))
}

fn which(name: &str) -> Option<PathBuf> {
    if let Ok(p) = std::env::var(format!("{}_PATH", name.to_ascii_uppercase())) {
        let pb = PathBuf::from(p);
        if pb.is_file() {
            return Some(pb);
        }
    }
    std::env::var_os("PATH").and_then(|paths| {
        std::env::split_paths(&paths).find_map(|dir| {
            let p = dir.join(name);
            p.is_file().then_some(p)
        })
    })
}

/// Strip Illumina `/1` `/2` mate suffixes for pairing across BAM and FASTQ.
pub fn normalize_qname(name: &str) -> &str {
    name.strip_suffix("/1")
        .or_else(|| name.strip_suffix("/2"))
        .unwrap_or(name)
}

/// Second-pass FASTQ: keep unmapped mates and reads with a seed-length clip.
pub fn keep_for_second_pass(unmapped: bool, mate_unmapped: bool, max_softclip: usize) -> bool {
    unmapped || mate_unmapped || max_softclip >= MIN_CLIP_FOR_SEED
}

/// Keep a FASTQ record if its qname is in `keep`, or never appeared in the BAM
/// (`seen`) — the latter are fully unmapped pairs dropped by `--no-unal`.
pub fn filter_fastq_file(
    src: &Path,
    dest: &Path,
    keep: &HashSet<String>,
    seen: &HashSet<String>,
) -> Result<usize> {
    if src
        .file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.ends_with(".gz"))
    {
        return Err(EvidenceError::Alignment(format!(
            "gzip FASTQ not supported for second-pass filter: {}",
            src.display()
        )));
    }
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let reader = BufReader::new(File::open(src)?);
    let mut lines = reader.lines();
    let mut out = File::create(dest)?;
    let mut n = 0usize;
    loop {
        let Some(header) = lines.next() else {
            break;
        };
        let header = header?;
        let seq = lines.next().ok_or_else(|| {
            EvidenceError::Alignment(format!("truncated FASTQ {}", src.display()))
        })??;
        let plus = lines.next().ok_or_else(|| {
            EvidenceError::Alignment(format!("truncated FASTQ {}", src.display()))
        })??;
        let qual = lines.next().ok_or_else(|| {
            EvidenceError::Alignment(format!("truncated FASTQ {}", src.display()))
        })??;
        let raw = header.strip_prefix('@').unwrap_or(header.as_str());
        let token = raw.split_whitespace().next().unwrap_or(raw);
        let q = normalize_qname(token);
        if keep.contains(q) || !seen.contains(q) {
            writeln!(out, "{header}")?;
            writeln!(out, "{seq}")?;
            writeln!(out, "{plus}")?;
            writeln!(out, "{qual}")?;
            n += 1;
        }
    }
    Ok(n)
}

/// Filter every file in `reads` into `dest_dir/pass_{i}.fastq`.
pub fn filter_fastq_input(
    reads: &FastqInput,
    keep: &HashSet<String>,
    seen: &HashSet<String>,
    dest_dir: &Path,
) -> Result<(FastqInput, usize)> {
    std::fs::create_dir_all(dest_dir)?;
    let mut files = Vec::new();
    let mut total = 0usize;
    for (i, src) in reads.files.iter().enumerate() {
        let dest = dest_dir.join(format!("pass_{i}.fastq"));
        total += filter_fastq_file(src, &dest, keep, seen)?;
        files.push(dest);
    }
    Ok((FastqInput { files }, total))
}

/// Bowtie2 `--un-conc prefix.fastq` writes `prefix.1.fastq` / `prefix.2.fastq`.
pub(crate) fn unconc_mate_paths(prefix: &Path) -> (PathBuf, PathBuf) {
    let file = prefix.file_name().and_then(|s| s.to_str()).unwrap_or("");
    let parent = prefix.parent().unwrap_or(Path::new(""));
    if let Some((stem, ext)) = file.rsplit_once('.') {
        (
            parent.join(format!("{stem}.1.{ext}")),
            parent.join(format!("{stem}.2.{ext}")),
        )
    } else {
        (
            parent.join(format!("{file}.1")),
            parent.join(format!("{file}.2")),
        )
    }
}

pub(crate) fn bowtie2_align_args(
    pass: AlignPass,
    threads: usize,
    index: &Path,
    sam: &Path,
    unmatched_conc: Option<&Path>,
    unmatched_se: Option<&Path>,
) -> Vec<String> {
    let mut args = vec![
        "--local".to_string(),
        "--no-unal".to_string(),
        "-p".to_string(),
        threads.max(1).to_string(),
        "-x".to_string(),
        index.to_string_lossy().into_owned(),
        "-S".to_string(),
        sam.to_string_lossy().into_owned(),
    ];
    match pass {
        AlignPass::PrimaryStage1 => {
            if let Some(p) = unmatched_conc {
                args.push("--un-conc".to_string());
                args.push(p.to_string_lossy().into_owned());
            }
            if let Some(p) = unmatched_se {
                args.push("--un".to_string());
                args.push(p.to_string_lossy().into_owned());
            }
        }
        AlignPass::PrimaryStage2 => {
            push_sensitive_bowtie2_opts(&mut args, UNMATCHED_SEED_LEN, "L,6,0.2", 200);
        }
        AlignPass::Junction => {
            push_sensitive_bowtie2_opts(&mut args, JUNCTION_SEED_LEN, "L,1,0.70", 2000);
        }
    }
    args
}

fn push_sensitive_bowtie2_opts(args: &mut Vec<String>, seed_len: usize, score_min: &str, k: usize) {
    args.extend([
        "--ma".to_string(),
        "1".to_string(),
        "--mp".to_string(),
        "3".to_string(),
        "--np".to_string(),
        "0".to_string(),
        "--rdg".to_string(),
        "2,3".to_string(),
        "--rfg".to_string(),
        "2,3".to_string(),
        "--ignore-quals".to_string(),
        "-L".to_string(),
        seed_len.to_string(),
        "-i".to_string(),
        "S,1,0.25".to_string(),
        "--score-min".to_string(),
        score_min.to_string(),
        "-k".to_string(),
        k.to_string(),
    ]);
}

fn append_fastq_args(args: &mut Vec<String>, reads: &FastqInput) -> Result<()> {
    match reads.files.as_slice() {
        [se] => {
            args.push("-U".to_string());
            args.push(se.to_string_lossy().into_owned());
        }
        [r1, r2] => {
            args.push("-1".to_string());
            args.push(r1.to_string_lossy().into_owned());
            args.push("-2".to_string());
            args.push(r2.to_string_lossy().into_owned());
        }
        files if files.len() >= 2 && files.len() % 2 == 0 => {
            let mut r1: Vec<String> = Vec::new();
            let mut r2: Vec<String> = Vec::new();
            for chunk in files.chunks(2) {
                r1.push(chunk[0].to_string_lossy().into_owned());
                r2.push(chunk[1].to_string_lossy().into_owned());
            }
            args.push("-1".to_string());
            args.push(r1.join(","));
            args.push("-2".to_string());
            args.push(r2.join(","));
        }
        extra => {
            return Err(EvidenceError::Alignment(format!(
                "expected 1 FASTQ (SE) or an even number of FASTQ files (PE pairs), got {}",
                extra.len()
            )));
        }
    }
    Ok(())
}

fn fastq_nonempty(path: &Path) -> bool {
    std::fs::metadata(path)
        .map(|m| m.len() > 0)
        .unwrap_or(false)
}

/// Append alignment lines from `extra` onto `base`, skipping `@` headers and
/// records whose (qname, mate) already appears in `base`. Stage 2 `--un-conc`
/// can contain pairs that already aligned discordantly in Stage 1; PE mates
/// share a qname so the key must include the mate bits, not qname alone.
#[cfg(test)]
pub(crate) fn append_sam_alignments(base: &Path, extra: &Path) -> Result<()> {
    let mut seen = sam_alignment_keys(base)?;
    let extra_f = File::open(extra)?;
    let mut out = std::fs::OpenOptions::new().append(true).open(base)?;
    for line in BufReader::new(extra_f).lines() {
        let line = line?;
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        if let Some(key) = sam_record_key(&line) {
            if !seen.insert(key) {
                continue;
            }
        }
        writeln!(out, "{line}")?;
    }
    Ok(())
}

#[cfg(test)]
fn sam_record_key(line: &str) -> Option<(String, u8)> {
    let mut fields = line.split('\t');
    let q = normalize_qname(fields.next()?);
    let flag: u16 = fields.next()?.parse().ok()?;
    let mate = if flag & 0x80 != 0 {
        2
    } else if flag & 0x40 != 0 {
        1
    } else {
        0
    };
    Some((q.to_string(), mate))
}

#[cfg(test)]
fn sam_alignment_keys(path: &Path) -> Result<HashSet<(String, u8)>> {
    let mut set = HashSet::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        if let Some(key) = sam_record_key(&line) {
            set.insert(key);
        }
    }
    Ok(set)
}

fn run_bowtie2(bin: &Path, args: &[String]) -> Result<()> {
    let out = Command::new(bin).args(args).output()?;
    if !out.status.success() {
        return Err(EvidenceError::Bowtie2 {
            status: out.status,
            stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
        });
    }
    Ok(())
}

/// `bowtie2-build` + `bowtie2 --local` (kind-specific seeds), then noodles SAM→sorted BAM.
pub fn align_to_bam(
    ref_fa: &Path,
    reads: &FastqInput,
    out_bam: &Path,
    threads: usize,
    work: &Path,
    kind: AlignKind,
) -> Result<()> {
    let bowtie2 = which_or_err("bowtie2")?;
    let build = which("bowtie2-build").unwrap_or_else(|| {
        bowtie2
            .parent()
            .map(|d| d.join("bowtie2-build"))
            .unwrap_or_else(|| PathBuf::from("bowtie2-build"))
    });
    if !build.is_file() && which("bowtie2-build").is_none() {
        return Err(EvidenceError::MissingBinary("bowtie2-build"));
    }

    std::fs::create_dir_all(work)?;
    let prefix = work.join("bt2ref");
    let status = Command::new(&build)
        .args([
            "--threads",
            &threads.max(1).to_string(),
            ref_fa.to_str().unwrap_or(""),
            prefix.to_str().unwrap_or(""),
        ])
        .output()?;
    if !status.status.success() {
        return Err(EvidenceError::Bowtie2Build {
            status: status.status,
            stderr: String::from_utf8_lossy(&status.stderr).into_owned(),
        });
    }

    let is_pe = reads.files.len() != 1;
    let sam_path = work.join("aligned.sam");
    let unconc = work.join("unconc.fastq");
    let unpaired = work.join("unpaired.fastq");
    let pass = match kind {
        AlignKind::Primary => AlignPass::PrimaryStage1,
        AlignKind::Junction => AlignPass::Junction,
    };
    let (unmatched_conc, unmatched_se) = match kind {
        AlignKind::Primary if is_pe => (Some(unconc.as_path()), Some(unpaired.as_path())),
        AlignKind::Primary => (None, Some(unpaired.as_path())),
        AlignKind::Junction => (None, None),
    };
    let mut args = bowtie2_align_args(
        pass,
        threads,
        &prefix,
        &sam_path,
        unmatched_conc,
        unmatched_se,
    );
    append_fastq_args(&mut args, reads)?;
    run_bowtie2(&bowtie2, &args)?;

    if kind == AlignKind::Primary {
        if is_pe {
            let (r1, r2) = unconc_mate_paths(&unconc);
            if fastq_nonempty(&r1) {
                eprintln!(
                    "prokdiff: primary stage-2 unmatched R1 (-L {UNMATCHED_SEED_LEN}, --score-min L,6,0.2)"
                );
                let sam_r1 = work.join("stage2_r1.sam");
                let mut s2 = bowtie2_align_args(
                    AlignPass::PrimaryStage2,
                    threads,
                    &prefix,
                    &sam_r1,
                    None,
                    None,
                );
                s2.push("-U".to_string());
                s2.push(r1.to_string_lossy().into_owned());
                run_bowtie2(&bowtie2, &s2)?;
            }
            if fastq_nonempty(&r2) {
                eprintln!(
                    "prokdiff: primary stage-2 unmatched R2 (-L {UNMATCHED_SEED_LEN}, --score-min L,6,0.2)"
                );
                let sam_r2 = work.join("stage2_r2.sam");
                let mut s2 = bowtie2_align_args(
                    AlignPass::PrimaryStage2,
                    threads,
                    &prefix,
                    &sam_r2,
                    None,
                    None,
                );
                s2.push("-U".to_string());
                s2.push(r2.to_string_lossy().into_owned());
                run_bowtie2(&bowtie2, &s2)?;
            }
        }
        if fastq_nonempty(&unpaired) {
            eprintln!(
                "prokdiff: primary stage-2 unmatched SE as single reads (-L {UNMATCHED_SEED_LEN}, --score-min L,6,0.2)"
            );
            let sam2 = work.join("stage2_se.sam");
            let mut s2 = bowtie2_align_args(
                AlignPass::PrimaryStage2,
                threads,
                &prefix,
                &sam2,
                None,
                None,
            );
            s2.push("-U".to_string());
            s2.push(unpaired.to_string_lossy().into_owned());
            run_bowtie2(&bowtie2, &s2)?;
            // stage2_se.sam is preserved in work for split candidate extraction.
        }
    }

    sam_to_sorted_bam(&sam_path, out_bam)?;
    Ok(())
}

pub(crate) fn sam_to_sorted_bam(sam_path: &Path, bam_path: &Path) -> Result<()> {
    let mut reader = File::open(sam_path)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;
    let header = reader
        .read_header()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    let mut records: Vec<sam::Record> = Vec::new();
    for rec in reader.records() {
        records.push(rec.map_err(|e| EvidenceError::Alignment(e.to_string()))?);
    }
    records.sort_by(|a, b| {
        let aid = a
            .reference_sequence_id(&header)
            .and_then(std::result::Result::ok);
        let bid = b
            .reference_sequence_id(&header)
            .and_then(std::result::Result::ok);
        let apos = a
            .alignment_start()
            .and_then(std::result::Result::ok)
            .map(|p| p.get())
            .unwrap_or(0);
        let bpos = b
            .alignment_start()
            .and_then(std::result::Result::ok)
            .map(|p| p.get())
            .unwrap_or(0);
        aid.cmp(&bid).then(apos.cmp(&bpos))
    });

    if let Some(parent) = bam_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut writer = File::create(bam_path).map(bam::io::Writer::new)?;
    writer
        .write_header(&header)
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    for rec in &records {
        writer
            .write_alignment_record(&header, rec)
            .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    }
    writer
        .try_finish()
        .map_err(|e| EvidenceError::Alignment(e.to_string()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn pe_input_keeps_r1_then_r2_order() {
        let fq = FastqInput::pe("a_R1.fastq", "a_R2.fastq");
        assert_eq!(fq.files.len(), 2);
        assert!(fq.files[0].ends_with("a_R1.fastq"));
        assert!(fq.files[1].ends_with("a_R2.fastq"));
    }

    #[test]
    fn normalize_qname_strips_mate_suffix() {
        assert_eq!(normalize_qname("SRR1.1"), "SRR1.1");
        assert_eq!(normalize_qname("SRR1.1/1"), "SRR1.1");
        assert_eq!(normalize_qname("SRR1.1/2"), "SRR1.1");
    }

    #[test]
    fn keep_for_second_pass_is_clip_or_unmapped() {
        assert!(!keep_for_second_pass(false, false, 0));
        assert!(!keep_for_second_pass(false, false, 5));
        assert!(keep_for_second_pass(false, false, 6));
        assert!(keep_for_second_pass(true, false, 0));
        assert!(keep_for_second_pass(false, true, 0));
    }

    #[test]
    fn primary_stage1_writes_unmatched_and_has_no_short_seed() {
        let args = bowtie2_align_args(
            AlignPass::PrimaryStage1,
            8,
            Path::new("idx"),
            Path::new("out.sam"),
            Some(Path::new("unconc.fastq")),
            None,
        );
        assert!(args.iter().any(|a| a == "--un-conc"));
        assert!(!args.windows(2).any(|w| w[0] == "-L"));
    }

    #[test]
    fn primary_stage2_uses_unmatched_seed_len() {
        let args = bowtie2_align_args(
            AlignPass::PrimaryStage2,
            8,
            Path::new("idx"),
            Path::new("out.sam"),
            None,
            None,
        );
        let want = UNMATCHED_SEED_LEN.to_string();
        assert!(
            args.windows(2).any(|w| w[0] == "-L" && w[1] == want),
            "stage-2 unmatched vs genome must use a short seed, got {args:?}"
        );
        assert!(
            args.windows(2)
                .any(|w| w[0] == "--score-min" && w[1] == "L,6,0.2"),
            "stage-2 must use sensitive score-min L,6,0.2; got {args:?}"
        );
        assert!(
            args.windows(2).any(|w| w[0] == "-k" && w[1] == "200"),
            "stage-2 must allow multi-hits (-k 200); got {args:?}"
        );
    }

    #[test]
    fn junction_second_pass_uses_seed_len_10() {
        let args = bowtie2_align_args(
            AlignPass::Junction,
            8,
            Path::new("idx"),
            Path::new("out.sam"),
            None,
            None,
        );
        assert!(
            args.windows(2).any(|w| w[0] == "-L" && w[1] == "10"),
            "36 bp reads cannot seed ~18 bp flanks at bowtie2 default -L 20; got {args:?}"
        );
        assert!(
            args.windows(2)
                .any(|w| w[0] == "--score-min" && w[1] == "L,1,0.70"),
            "junction second pass must use score-min L,1,0.70; got {args:?}"
        );
        assert!(
            args.windows(2).any(|w| w[0] == "-k" && w[1] == "2000"),
            "junction second pass must allow multi-hits (-k 2000); got {args:?}"
        );
        assert!(!args.iter().any(|a| a == "--un-conc" || a == "--un"));
    }

    #[test]
    fn unconc_paths_insert_mate_before_fastq_suffix() {
        let (r1, r2) = unconc_mate_paths(Path::new("/tmp/unconc.fastq"));
        assert!(r1.ends_with("unconc.1.fastq"));
        assert!(r2.ends_with("unconc.2.fastq"));
    }

    #[test]
    fn append_sam_skips_extra_headers() {
        let dir = std::env::temp_dir().join("prokdiff-sam-merge");
        let _ = std::fs::create_dir_all(&dir);
        let base = dir.join("a.sam");
        let extra = dir.join("b.sam");
        std::fs::write(
            &base,
            "@HD\tVN:1.0\nr1\t0\tref\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();
        std::fs::write(
            &extra,
            "@HD\tVN:1.0\n@SQ\tSN:ref\tLN:4\nr1\t0\tref\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\nr2\t0\tref\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();
        append_sam_alignments(&base, &extra).unwrap();
        let out = std::fs::read_to_string(&base).unwrap();
        assert_eq!(out.matches("@HD").count(), 1);
        assert_eq!(out.matches("r1\t").count(), 1);
        assert_eq!(out.matches("r2\t").count(), 1);
    }

    #[test]
    fn append_sam_keeps_both_mates_of_a_new_pair() {
        let dir = std::env::temp_dir().join("prokdiff-sam-merge-pe");
        let _ = std::fs::create_dir_all(&dir);
        let base = dir.join("a.sam");
        let extra = dir.join("b.sam");
        std::fs::write(&base, "@HD\tVN:1.0\n").unwrap();
        std::fs::write(
            &extra,
            "pair\t99\tref\t1\t255\t4M\t=\t10\t0\tACGT\tIIII\npair\t147\tref\t10\t255\t4M\t=\t1\t0\tACGT\tIIII\n",
        )
        .unwrap();
        append_sam_alignments(&base, &extra).unwrap();
        let out = std::fs::read_to_string(&base).unwrap();
        assert_eq!(out.matches("pair\t").count(), 2);
    }

    #[test]
    fn filter_fastq_keeps_clipped_and_unseen_drops_mapped_only() {
        let dir = std::env::temp_dir().join("prokdiff-fq-filter");
        let _ = std::fs::create_dir_all(&dir);
        let src = dir.join("in.fastq");
        std::fs::write(
            &src,
            "@clip1\nACGT\n+\nIIII\n@mapped_only\nACGT\n+\nIIII\n@unmapped\nTTTT\n+\nIIII\n",
        )
        .unwrap();
        let dest = dir.join("out.fastq");
        let keep = HashSet::from(["clip1".to_string()]);
        let seen = HashSet::from(["clip1".to_string(), "mapped_only".to_string()]);
        let n = filter_fastq_file(&src, &dest, &keep, &seen).unwrap();
        assert_eq!(n, 2);
        let out = std::fs::read_to_string(&dest).unwrap();
        assert!(out.contains("@clip1"));
        assert!(out.contains("@unmapped"));
        assert!(!out.contains("@mapped_only"));
    }
}
