#![deny(unsafe_code)]

pub mod associate;
pub mod io;
pub mod model;
pub mod normalize;
pub mod oracle;
pub mod rust_search;
pub mod scoring;

pub use associate::link_mutations_to_sites;
pub use io::{write_mutation_offtarget_links_tsv, write_offtarget_sites_tsv};
pub use model::{
    BulgeType, MutationOffTargetLink, NucleaseProfile, OffTargetSite, PamSide, Strand,
};
pub use oracle::{parse_cas_offinder, parse_crispritz, parse_flashfry};
pub use rust_search::{
    scan_contig, scan_contig_bulge, scan_genome, scan_genome_bulge, BulgeSearchOptions,
};
pub use scoring::{calculate_cfd_score, calculate_hsu_score};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_cas_offinder_output() {
        let sample = "\
GAGTCCGAGCAGAAGAAGAA	chr1	100	GAGTCCGAGCAGAAGAAGAATGG	+	0
GAGTCCGAGCAGAAGAAGAA	chr1	250	GAGTCCGAGCAGAAGAACAAGGG	-	1
";
        let sites = parse_cas_offinder(sample).expect("should parse Cas-OFFinder output");
        assert_eq!(sites.len(), 2);

        // Site 1: 0-based 100 -> 1-based 101, len 23 -> end 123
        assert_eq!(sites[0].site_id, "cas_offinder_1");
        assert_eq!(sites[0].seq_id, "chr1");
        assert_eq!(sites[0].start, 101);
        assert_eq!(sites[0].end, 123);
        assert_eq!(sites[0].strand, Strand::Plus);
        assert_eq!(sites[0].mismatches, 0);
        assert_eq!(sites[0].pam, "TGG");
        assert_eq!(sites[0].search_backend, "cas-offinder");

        // Site 2: 0-based 250 -> 1-based 251, len 23 -> end 273
        assert_eq!(sites[1].site_id, "cas_offinder_2");
        assert_eq!(sites[1].start, 251);
        assert_eq!(sites[1].end, 273);
        assert_eq!(sites[1].strand, Strand::Minus);
        assert_eq!(sites[1].mismatches, 1);
        assert_eq!(sites[1].pam, "GGG");
    }

    #[test]
    fn parses_flashfry_output() {
        let sample = "\
contig	start	stop	strand	target	mismatches	cfd	hsu
NC_000913.3	500	523	+	GAGTCCGAGCAGAAGAAGAATGG	0	1.0000	99.5
NC_000913.3	1200	1223	-	GAGTCCGAGCAGAAGAACAACGG	1	0.7500	82.1
";
        let sites = parse_flashfry(sample).expect("should parse FlashFry output");
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].seq_id, "NC_000913.3");
        assert_eq!(sites[0].start, 500);
        assert_eq!(sites[0].end, 523);
        assert_eq!(sites[0].strand, Strand::Plus);
        assert_eq!(sites[0].mismatches, 0);
        assert_eq!(sites[0].cfd_score, Some(1.0));
        assert_eq!(sites[0].hsu_score, Some(99.5));

        assert_eq!(sites[1].strand, Strand::Minus);
        assert_eq!(sites[1].mismatches, 1);
        assert_eq!(sites[1].cfd_score, Some(0.75));
    }

    #[test]
    fn parses_crispritz_output() {
        let sample = "\
chr	position	stop	target	pam	mismatches	bulge_type	bulge_size	strand	cfd
chr1	1000	1023	GAGTCCGAGCAGAAGAAGAATGG	NGG	0	none	0	+	1.0
chr1	2000	2024	GAGTCCGAGCAGAAGAACAATGG	NGG	1	dna	1	-	0.55
";
        let sites = parse_crispritz(sample).expect("should parse CRISPRitz output");
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].start, 1000);
        assert_eq!(sites[0].end, 1023);
        assert_eq!(sites[0].strand, Strand::Plus);
        assert_eq!(sites[0].bulge_type, BulgeType::None);

        assert_eq!(sites[1].strand, Strand::Minus);
        assert_eq!(sites[1].bulge_type, BulgeType::Dna);
        assert_eq!(sites[1].bulge_size, 1);
        assert_eq!(sites[1].cfd_score, Some(0.55));
    }

    #[test]
    fn rust_exact_search_finds_plus_and_minus_cas9_sites() {
        let spacer = "GAGTCCGAGCAGAAGAAGAA"; // 20 bp
        let profile = NucleaseProfile::spcas9(); // 3' NGG

        // Construct synthetic reference:
        // Index 10: spacer + "CGG" (Plus match at start 11, end 33, 0 mm)
        // Index 50: revcomp("AGG") + revcomp(spacer) (Minus match at start 51, end 73, 0 mm)
        let mut ref_seq = vec![b'A'; 100];
        let plus_target = format!("{}CGG", spacer);
        ref_seq[10..10 + 23].copy_from_slice(plus_target.as_bytes());

        let pam_agg_rc = b"CCT"; // revcomp of AGG
        let spacer_rc = normalize::revcomp_dna(spacer.as_bytes());
        ref_seq[50..53].copy_from_slice(pam_agg_rc);
        ref_seq[53..53 + 20].copy_from_slice(&spacer_rc);

        let sites = scan_contig("synth_chr", &ref_seq, spacer, &profile, 0, 0);
        assert_eq!(sites.len(), 2);

        // Plus site
        let s_plus = sites.iter().find(|s| s.strand == Strand::Plus).unwrap();
        assert_eq!(s_plus.start, 11);
        assert_eq!(s_plus.end, 33);
        assert_eq!(s_plus.mismatches, 0);
        assert_eq!(s_plus.pam, "CGG");

        // Minus site
        let s_minus = sites.iter().find(|s| s.strand == Strand::Minus).unwrap();
        assert_eq!(s_minus.start, 51);
        assert_eq!(s_minus.end, 73);
        assert_eq!(s_minus.mismatches, 0);
        assert_eq!(s_minus.pam, "AGG");
    }

    #[test]
    fn rust_exact_search_finds_cas12a_sites_with_tttv_pam() {
        let spacer = "ACCAGTCAGTCAGTCAGTCA"; // 20 bp
        let profile = NucleaseProfile::cas12a(); // 5' TTTV (V = A, C, G)

        // Cas12a 5' PAM: [TTTG][spacer]
        let mut ref_seq = vec![b'C'; 100];
        let plus_target = format!("TTTG{}", spacer);
        ref_seq[20..20 + 24].copy_from_slice(plus_target.as_bytes());

        let sites = scan_contig("synth_cas12a", &ref_seq, spacer, &profile, 0, 0);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].strand, Strand::Plus);
        assert_eq!(sites[0].start, 21);
        assert_eq!(sites[0].end, 44);
        assert_eq!(sites[0].pam, "TTTG");
        assert_eq!(sites[0].mismatches, 0);
    }

    #[test]
    fn rust_exact_matches_cas_offinder_format_output() {
        let spacer = "GAGTCCGAGCAGAAGAAGAA";
        let profile = NucleaseProfile::spcas9();

        let mut ref_seq = vec![b'A'; 60];
        let plus_target = format!("{}CGG", spacer);
        ref_seq[10..10 + 23].copy_from_slice(plus_target.as_bytes());

        let rust_sites = scan_contig("chr1", &ref_seq, spacer, &profile, 2, 0);

        // Simulated Cas-OFFinder output line for this sequence:
        let cas_offinder_tsv = format!("{spacer}\tchr1\t10\t{plus_target}\t+\t0\n");
        let oracle_sites = parse_cas_offinder(&cas_offinder_tsv).unwrap();

        assert_eq!(rust_sites.len(), 1);
        assert_eq!(oracle_sites.len(), 1);

        assert_eq!(rust_sites[0].seq_id, oracle_sites[0].seq_id);
        assert_eq!(rust_sites[0].start, oracle_sites[0].start);
        assert_eq!(rust_sites[0].end, oracle_sites[0].end);
        assert_eq!(rust_sites[0].strand, oracle_sites[0].strand);
        assert_eq!(rust_sites[0].mismatches, oracle_sites[0].mismatches);
    }
}
