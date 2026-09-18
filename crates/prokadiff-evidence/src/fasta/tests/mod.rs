use super::*;
use std::io::Write;

fn write_fixture(path: &Path, contents: &str) {
    let mut file = File::create(path).unwrap();
    file.write_all(contents.as_bytes()).unwrap();
}

#[test]
fn reads_simple_fasta() {
    let dir = std::env::temp_dir().join("prokdiff-fasta-test");
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("t.fa");
    {
        let mut f = File::create(&path).unwrap();
        writeln!(f, ">chr comment\nACGT\nacgt").unwrap();
    }
    let recs = read_reference(&path).unwrap();
    assert_eq!(recs.len(), 1);
    assert_eq!(recs[0].name, "chr");
    assert_eq!(recs[0].seq, b"ACGTACGT");
}

#[test]
fn reads_genbank_origin() {
    let dir = std::env::temp_dir().join("prokdiff-fasta-test");
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("t.gbk");
    {
        let mut f = File::create(&path).unwrap();
        writeln!(f, "LOCUS       syn  8 bp\nORIGIN\n        1 acgtacgt\n//").unwrap();
    }
    let recs = read_reference(&path).unwrap();
    assert_eq!(recs[0].name, "syn");
    assert_eq!(recs[0].seq, b"ACGTACGT");
}

#[test]
fn reads_genbank_origin_prefers_version_over_locus_accession() {
    // RW-005: LOCUS carries the bare accession (no version suffix), while
    // VERSION carries the versioned accession that FASTA-derived seq_ids and
    // `--intended` TSVs actually use. The contig name must come from VERSION.
    let dir = std::env::temp_dir().join("prokdiff-fasta-test");
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("versioned.gbk");
    {
        let mut f = File::create(&path).unwrap();
        writeln!(
            f,
            "LOCUS       foo  8 bp\nVERSION     foo.2\nORIGIN\n        1 acgtacgt\n//"
        )
        .unwrap();
    }
    let recs = read_reference(&path).unwrap();
    assert_eq!(recs[0].name, "foo.2");
    assert_eq!(recs[0].seq, b"ACGTACGT");
}

#[test]
fn read_reference_rejects_duplicate_ids_in_one_multirecord_fasta() {
    // Given: one multi-record FASTA whose records share an identifier.
    let dir = std::env::temp_dir().join(format!(
        "prokadiff-fasta-duplicate-single-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let source = dir.join("duplicate.fa");
    write_fixture(&source, ">chr\nACGT\n>chr second copy\nTGCA\n");

    // When: the reference boundary reads the single input path.
    let error = read_reference(&source).expect_err("duplicate ID must be rejected");

    // Then: the typed error identifies the repeated identifier and both source locations.
    match error {
        EvidenceError::DuplicateReferenceId {
            id,
            first_path,
            duplicate_path,
        } => {
            assert_eq!(id, "chr");
            assert_eq!(first_path, source);
            assert_eq!(duplicate_path, source);
        }
        other => panic!("expected duplicate reference ID error, got {other}"),
    }
}

#[test]
fn read_reference_rejects_duplicate_ids_in_one_multirecord_genbank() {
    // Given: one multi-record GenBank input whose records share an identifier.
    let dir = std::env::temp_dir().join(format!(
        "prokadiff-genbank-duplicate-single-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let source = dir.join("duplicate.gbk");
    write_fixture(
        &source,
        "LOCUS       chr  4 bp\nORIGIN\n        1 acgt\n//\nLOCUS       chr  4 bp\nORIGIN\n        1 tgca\n//\n",
    );

    // When: the reference boundary reads the single input path.
    let error = read_reference(&source).expect_err("duplicate ID must be rejected");

    // Then: the same typed duplicate-ID error is returned.
    assert!(matches!(error, EvidenceError::DuplicateReferenceId { .. }));
}

#[test]
fn write_combined_fasta_rejects_duplicate_ids_across_paths() {
    // Given: two reference files that define the same contig identifier.
    let dir = std::env::temp_dir().join(format!(
        "prokadiff-fasta-duplicate-cross-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let first = dir.join("first.fa");
    let second = dir.join("second.fa");
    let dest = dir.join("combined.fa");
    write_fixture(&first, ">chr\nACGT\n");
    write_fixture(&second, ">chr\nTGCA\n");

    // When: the files are combined for an aligner reference.
    let error = write_combined_fasta(&[&first, &second], &dest)
        .expect_err("duplicate ID must be rejected before writing");

    // Then: the typed error identifies both conflicting paths.
    match error {
        EvidenceError::DuplicateReferenceId {
            id,
            first_path,
            duplicate_path,
        } => {
            assert_eq!(id, "chr");
            assert_eq!(first_path, first);
            assert_eq!(duplicate_path, second);
        }
        other => panic!("expected duplicate reference ID error, got {other}"),
    }
    assert!(!dest.exists());
}

#[test]
fn write_combined_fasta_accepts_unique_ids_across_paths() {
    // Given: two reference files whose contig identifiers are distinct.
    let dir = std::env::temp_dir().join(format!("prokadiff-fasta-unique-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let first = dir.join("first.fa");
    let second = dir.join("second.fa");
    let dest = dir.join("combined.fa");
    write_fixture(&first, ">chromosome\nACGT\n");
    write_fixture(&second, ">plasmid\nTGCA\n");

    // When: the files are combined for an aligner reference.
    write_combined_fasta(&[&first, &second], &dest).unwrap();

    // Then: both distinct contigs are present in the written FASTA.
    let combined = std::fs::read_to_string(&dest).unwrap();
    assert!(combined.contains(">chromosome\nACGT\n"));
    assert!(combined.contains(">plasmid\nTGCA\n"));
}

#[test]
fn parses_genbank_repeat_features() {
    let dir = std::env::temp_dir().join("prokdiff-fasta-test");
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("rep.gbk");
    {
        let mut f = File::create(&path).unwrap();
        writeln!(
            f,
            "LOCUS       chr1  1000 bp\nFEATURES             Location/Qualifiers\n     repeat_region   100..200\n                     /mobile_element=\"insertion sequence:IS150\"\n     repeat_region   complement(500..600)\n                     /mobile_element=\"IS186\"\nORIGIN\n        1 aaaa\n//"
        )
        .unwrap();
    }
    let reps = parse_genbank_repeats(&path).unwrap();
    assert_eq!(reps.len(), 2);
    assert_eq!(reps[0].seq_id, "chr1");
    assert_eq!(reps[0].start, 100);
    assert_eq!(reps[0].end, 200);
    assert_eq!(reps[0].strand, 1);
    assert_eq!(reps[0].name, "IS150");

    assert_eq!(reps[1].seq_id, "chr1");
    assert_eq!(reps[1].start, 500);
    assert_eq!(reps[1].end, 600);
    assert_eq!(reps[1].strand, -1);
    assert_eq!(reps[1].name, "IS186");
}

#[test]
fn parses_rel606_gbk_repeats() {
    let path = Path::new("../../testdata/layer2/clonal/Clonal_Sample/REL606.gbk");
    if !path.is_file() {
        return;
    }
    let reps = parse_genbank_repeats(path).unwrap();
    assert!(!reps.is_empty());
    let is150 = reps.iter().filter(|r| r.name == "IS150").count();
    let is186 = reps.iter().filter(|r| r.name == "IS186").count();
    assert!(is150 >= 5, "expected >=5 IS150 copies, found {}", is150);
    assert!(is186 >= 5, "expected >=5 IS186 copies, found {}", is186);
}

#[test]
fn parses_genbank_repeats_prefers_version_over_locus_accession() {
    // RW-005: repeat_region/mobile_element seq_id must also use VERSION.
    let dir = std::env::temp_dir().join("prokdiff-fasta-test");
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join("rep_versioned.gbk");
    {
        let mut f = File::create(&path).unwrap();
        writeln!(
            f,
            "LOCUS       chr1  1000 bp\nVERSION     chr1.3\nFEATURES             Location/Qualifiers\n     repeat_region   100..200\n                     /mobile_element=\"insertion sequence:IS150\"\nORIGIN\n        1 aaaa\n//"
        )
        .unwrap();
    }
    let reps = parse_genbank_repeats(&path).unwrap();
    assert_eq!(reps.len(), 1);
    assert_eq!(reps[0].seq_id, "chr1.3");
}

#[test]
fn test_parses_genbank_features() {
    use crate::fasta::parse_genbank_features;
    let dir = std::env::temp_dir().join(format!("test_gbk_feat_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("feat.gbk");
    {
        let mut f = File::create(&path).unwrap();
        write!(
            f,
            "LOCUS       chr1  1000 bp\nFEATURES             Location/Qualifiers\n     CDS             100..300\n                     /gene=\"dnaA\"\n                     /locus_tag=\"b0001\"\n                     /product=\"replication initiator\"\n     tRNA            complement(400..480)\n                     /gene=\"tRNA-Ala\"\n                     /locus_tag=\"b0002\"\nORIGIN\n        1 aaaa\n//\n"
        )
        .unwrap();
    }
    let feats = parse_genbank_features(&path).unwrap();
    assert_eq!(feats.len(), 2);
    assert_eq!(feats[0].seq_id, "chr1");
    assert_eq!(feats[0].start, 100);
    assert_eq!(feats[0].end, 300);
    assert_eq!(feats[0].strand, 1);
    assert_eq!(feats[0].feature_type, "CDS");
    assert_eq!(feats[0].gene_name.as_deref(), Some("dnaA"));
    assert_eq!(feats[0].locus_tag.as_deref(), Some("b0001"));
    assert_eq!(feats[0].product.as_deref(), Some("replication initiator"));

    assert_eq!(feats[1].seq_id, "chr1");
    assert_eq!(feats[1].start, 400);
    assert_eq!(feats[1].end, 480);
    assert_eq!(feats[1].strand, -1);
    assert_eq!(feats[1].feature_type, "tRNA");
    assert_eq!(feats[1].gene_name.as_deref(), Some("tRNA-Ala"));

    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn parses_genbank_features_prefers_version_over_locus_accession() {
    // RW-005: CDS/tRNA/rRNA/gene seq_id must also use VERSION.
    use crate::fasta::parse_genbank_features;
    let dir = std::env::temp_dir().join(format!("test_gbk_feat_ver_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("feat_versioned.gbk");
    {
        let mut f = File::create(&path).unwrap();
        write!(
            f,
            "LOCUS       chr1  1000 bp\nVERSION     chr1.5\nFEATURES             Location/Qualifiers\n     CDS             100..300\n                     /gene=\"dnaA\"\nORIGIN\n        1 aaaa\n//\n"
        )
        .unwrap();
    }
    let feats = parse_genbank_features(&path).unwrap();
    assert_eq!(feats.len(), 1);
    assert_eq!(feats[0].seq_id, "chr1.5");

    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn test_parses_multirecord_genbank_and_complex_features() {
    let dir = std::env::temp_dir().join(format!("test_gbk_multi_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("multi.gbk");
    {
        let mut f = File::create(&path).unwrap();
        write!(
            f,
            r#"LOCUS       NC_chr       1000 bp    DNA     circular BCT 16-SEP-2026
FEATURES             Location/Qualifiers
     gene            <100..>400
                     /gene="lacZ"
                     /locus_tag="b0344"
     CDS             join(100..200,300..400)
                     /gene="lacZ"
                     /locus_tag="b0344"
                     /product="beta-galactosidase
                     subunit alpha"
     repeat_region   500..600
                     /mobile_element="insertion sequence:IS1"
ORIGIN
        1 aaaaaaaa
//
LOCUS       pPlasmid_1    500 bp    DNA     circular BCT 16-SEP-2026
FEATURES             Location/Qualifiers
     CDS             complement(50..150)
                     /gene="repA"
                     /locus_tag="p0001"
                     /product="plasmid replication protein"
     repeat_region   complement(200..300)
                     /mobile_element="IS186"
ORIGIN
        1 tttttttt
//
"#
        )
        .unwrap();
    }

    let feats = parse_genbank_features(&path).unwrap();
    assert_eq!(feats.len(), 3, "should parse 3 features across 2 records");

    // Feature 1: gene on chr
    assert_eq!(feats[0].seq_id, "NC_chr");
    assert_eq!(feats[0].start, 100);
    assert_eq!(feats[0].end, 400);
    assert_eq!(feats[0].feature_type, "gene");
    assert_eq!(feats[0].gene_name.as_deref(), Some("lacZ"));

    // Feature 2: joined CDS on chr with multiline product
    assert_eq!(feats[1].seq_id, "NC_chr");
    assert_eq!(feats[1].start, 100);
    assert_eq!(feats[1].end, 400);
    assert_eq!(feats[1].feature_type, "CDS");
    assert_eq!(
        feats[1].product.as_deref(),
        Some("beta-galactosidase subunit alpha")
    );

    // Feature 3: CDS on plasmid (proves multi-record past first ORIGIN works!)
    assert_eq!(feats[2].seq_id, "pPlasmid_1");
    assert_eq!(feats[2].start, 50);
    assert_eq!(feats[2].end, 150);
    assert_eq!(feats[2].strand, -1);
    assert_eq!(feats[2].gene_name.as_deref(), Some("repA"));

    // Check repeats across both records
    let reps = parse_genbank_repeats(&path).unwrap();
    assert_eq!(reps.len(), 2, "should parse 2 repeats across 2 records");
    assert_eq!(reps[0].seq_id, "NC_chr");
    assert_eq!(reps[0].name, "IS1");
    assert_eq!(reps[1].seq_id, "pPlasmid_1");
    assert_eq!(reps[1].name, "IS186");

    let _ = std::fs::remove_dir_all(dir);
}
