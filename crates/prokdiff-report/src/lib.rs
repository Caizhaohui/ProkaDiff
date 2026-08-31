//! unintended.tsv and a short run summary (stage 3 surface; stage 4 may extend).

#![deny(unsafe_code)]

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use prokdiff_classify::{ClassifiedMutation, ClassifyResult, MutationClass};
use prokdiff_gd::{GdEntry, GdKind};

/// Write classified post-subtract mutations. `hypothesis` column omitted when `include_hypothesis` is false.
pub fn write_unintended_tsv(
    path: impl AsRef<Path>,
    rows: &[ClassifiedMutation],
    editor: &str,
    include_hypothesis: bool,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write!(
        w,
        "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site"
    )?;
    if include_hypothesis {
        write!(w, "\thypothesis")?;
    }
    writeln!(w)?;
    for row in rows {
        let (seq_id, pos, end) = coords(&row.entry);
        let (ref_a, alt) = alleles(&row.entry);
        write!(
            w,
            "{seq_id}\t{pos}\t{end}\t{}\t{ref_a}\t{alt}\t{}\t{editor}\t{}\t{}\t{}",
            row.entry.kind.as_str(),
            row.class.as_str(),
            row.pam_profile.as_deref().unwrap_or(""),
            opt_u32(row.offtarget_mismatch),
            opt_u64(row.distance_to_site),
        )?;
        if include_hypothesis {
            write!(w, "\t{}", row.hypothesis.as_deref().unwrap_or(""))?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}

pub fn write_summary(
    path: impl AsRef<Path>,
    result: &ClassifyResult,
    intended_provided: bool,
    editor: &str,
) -> std::io::Result<()> {
    let mut n_s = 0usize;
    let mut n_n = 0usize;
    let mut n_c = 0usize;
    for r in &result.unintended {
        match r.class {
            MutationClass::Structural => n_s += 1,
            MutationClass::NearHomolog => n_n += 1,
            MutationClass::ScatteredSnv => n_c += 1,
        }
    }
    let intended_line = if intended_provided {
        format!("{}", result.intended_observed.len())
    } else {
        "NA".into()
    };
    let text = format!(
        "editor\t{editor}\n\
intended_provided\t{}\n\
intended_observed\t{intended_line}\n\
structural\t{n_s}\n\
near_homolog\t{n_n}\n\
scattered_snv\t{n_c}\n\
starter_vs_ref_mutations\t{}\n",
        if intended_provided { "yes" } else { "no" },
        result.starter_vs_ref,
    );
    std::fs::write(path, text)?;
    Ok(())
}

fn opt_u32(v: Option<u32>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}

fn opt_u64(v: Option<u64>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}

fn coords(e: &GdEntry) -> (String, u64, u64) {
    let seq = e.seq_id().unwrap_or("").to_string();
    let pos = e.position().unwrap_or(0);
    if e.kind == GdKind::Del {
        let size = e
            .fields
            .get(2)
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(1)
            .max(1);
        (seq, pos, pos.saturating_add(size - 1))
    } else {
        (seq, pos, pos)
    }
}

fn alleles(e: &GdEntry) -> (&str, &str) {
    match e.kind {
        GdKind::Snp | GdKind::Ins => (".", e.fields.get(2).map(String::as_str).unwrap_or(".")),
        _ => (".", "."),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use prokdiff_classify::{ClassifiedMutation, MutationClass};
    use prokdiff_gd::GdEntry;

    #[test]
    fn tsv_omits_hypothesis_column_when_disabled() {
        let dir = std::env::temp_dir().join("prokdiff-report-test");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("unintended.tsv");
        let row = ClassifiedMutation {
            entry: GdEntry::snp(1, "chr", 40, "C"),
            class: MutationClass::ScatteredSnv,
            pam_profile: None,
            offtarget_mismatch: None,
            distance_to_site: None,
            hypothesis: Some("sos_widney2014".into()),
        };
        write_unintended_tsv(&path, &[row], "cas9", false).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.starts_with(
            "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\n"
        ));
        assert!(!text.contains("hypothesis"));
        assert!(text.contains("scattered_snv"));
        assert!(text.contains("cas9"));
    }

    #[test]
    fn summary_marks_intended_na_when_omitted() {
        let dir = std::env::temp_dir().join("prokdiff-report-test");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("summary.txt");
        let result = ClassifyResult {
            unintended: vec![ClassifiedMutation {
                entry: GdEntry::snp(1, "chr", 40, "C"),
                class: MutationClass::NearHomolog,
                pam_profile: Some("NGG".into()),
                offtarget_mismatch: Some(0),
                distance_to_site: Some(7),
                hypothesis: None,
            }],
            intended_observed: vec![],
            starter_vs_ref: 50,
        };
        write_summary(&path, &result, false, "cas9").unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("intended_observed\tNA"));
        assert!(text.contains("near_homolog\t1"));
        assert!(text.contains("starter_vs_ref_mutations\t50"));
    }
}
