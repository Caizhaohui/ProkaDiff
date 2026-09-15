use prokadiff_offtarget::{
    scan_genome_bulge, write_offtarget_sites_tsv, BulgeSearchOptions, NucleaseProfile,
};
use std::path::{Path, PathBuf};

fn read_simple_fasta(path: &Path) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    let text = std::fs::read_to_string(path)?;
    let mut recs = Vec::new();
    let mut cur_name = String::new();
    let mut cur_seq = Vec::new();
    for line in text.lines() {
        let l = line.trim();
        if l.is_empty() {
            continue;
        }
        if let Some(name) = l.strip_prefix('>') {
            if !cur_name.is_empty() {
                recs.push((cur_name, cur_seq));
                cur_seq = Vec::new();
            }
            cur_name = name.split_whitespace().next().unwrap_or("").to_string();
        } else {
            cur_seq.extend_from_slice(l.as_bytes());
        }
    }
    if !cur_name.is_empty() {
        recs.push((cur_name, cur_seq));
    }
    Ok(recs)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut ref_paths = Vec::new();
    let mut spacer = String::new();
    let mut pam = String::from("NGG");
    let mut editor = String::from("cas9");
    let mut max_mismatches = 3u32;
    let mut max_dna_bulge = 0u32;
    let mut max_rna_bulge = 0u32;
    let mut output_path = PathBuf::from("offtarget_sites.tsv");

    let args: Vec<String> = std::env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--ref" => {
                ref_paths.push(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--spacer" => {
                spacer = args[i + 1].clone();
                i += 2;
            }
            "--pam" => {
                pam = args[i + 1].clone();
                i += 2;
            }
            "--editor" => {
                editor = args[i + 1].to_ascii_lowercase();
                i += 2;
            }
            "--max-mismatches" => {
                max_mismatches = args[i + 1].parse()?;
                i += 2;
            }
            "--max-dna-bulge" => {
                max_dna_bulge = args[i + 1].parse()?;
                i += 2;
            }
            "--max-rna-bulge" => {
                max_rna_bulge = args[i + 1].parse()?;
                i += 2;
            }
            "--output" => {
                output_path = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            _ => {
                eprintln!("Unknown argument: {}", args[i]);
                i += 1;
            }
        }
    }

    if ref_paths.is_empty() || spacer.is_empty() {
        eprintln!("Usage: offtarget_scan --ref <ref.fa> --spacer <guide> [options]");
        std::process::exit(1);
    }

    let mut refs = Vec::new();
    for p in &ref_paths {
        let recs = read_simple_fasta(p)?;
        refs.extend(recs);
    }

    let profile = match editor.as_str() {
        "cas12a" => NucleaseProfile::cas12a(),
        _ => {
            let mut prof = NucleaseProfile::spcas9();
            prof.pam_pattern = pam;
            prof
        }
    };

    let opts = BulgeSearchOptions {
        max_mismatches,
        max_dna_bulge,
        max_rna_bulge,
    };

    let sites = scan_genome_bulge(&refs, &spacer, &profile, opts);
    write_offtarget_sites_tsv(&sites, &output_path)?;

    println!(
        "Found {} off-target sites. Wrote {}",
        sites.len(),
        output_path.display()
    );
    Ok(())
}
